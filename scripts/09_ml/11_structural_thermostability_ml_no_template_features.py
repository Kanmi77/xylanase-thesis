#!/usr/bin/env python3

from pathlib import Path
import warnings
import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


warnings.filterwarnings("ignore")

BASE = Path.home() / "xylanase-thesis"

INFILE = BASE / "data/curated/xylanase_thesis_master_final_v5_activity_labels.csv"

OUTDIR = BASE / "results/ml"
REPORTDIR = BASE / "results/reports"

OUTDIR.mkdir(parents=True, exist_ok=True)
REPORTDIR.mkdir(parents=True, exist_ok=True)

PRED_OUT = OUTDIR / "structural_thermostability_ml_no_template_predictions.csv"
SUMMARY_OUT = OUTDIR / "structural_thermostability_ml_no_template_summary.csv"
FEATURE_OUT = OUTDIR / "structural_thermostability_ml_no_template_feature_importance.csv"
PERM_OUT = OUTDIR / "structural_thermostability_ml_no_template_permutation_importance.csv"
GROUP_OUT = OUTDIR / "structural_thermostability_ml_no_template_group_summary.csv"
AUDIT_OUT = REPORTDIR / "structural_thermostability_ml_no_template_audit.txt"


TARGET = "foldx_energy_per_residue"

# Template quality features deliberately excluded:
# - modeller_best_template_identity
# - modeller_best_template_coverage
#
# This sensitivity model tests whether structural thermostability prediction remains strong
# when model-template quality variables are removed.

NUMERIC_FEATURES = [
    # Sequence and physicochemical features
    "sequence_length",
    "molecular_weight",
    "aromaticity",
    "instability_index",
    "isoelectric_point",
    "gravy",
    "helix_fraction",
    "turn_fraction",
    "sheet_fraction",

    # Signal peptide proxy
    "nterm_hydrophobic_count_5_25",
    "nterm_positive_count_1_10",

    # Structural size feature
    "chain_length",

    # Structural stability-associated features
    "hbond_proxy_count",
    "salt_bridge_count",
    "disulfide_count",
    "sasa_total",
    "hbond_per_res",
    "salt_bridge_per_res",
    "disulfide_per_res",
    "sasa_per_res",
]

CATEGORICAL_FEATURES = [
    "gh_family",
    "organism_type",
    "structure_sources",
    "foldx_source_layer",
    "signal_peptide_proxy",
    "sequence_has_X",
    "sequence_has_nonstandard_residue",
    "is_short_sequence_lt100aa",
    "is_long_sequence_gt1200aa",
]


EXCLUDE_FROM_FEATURES = [
    # Direct target/leakage
    "foldx_wt_total_energy",
    "foldx_energy_per_residue",
    "predicted_foldx_energy_per_residue",
    "prediction_error",
    "has_ml_prediction",

    # Template-quality variables deliberately excluded in this sensitivity model
    "modeller_best_template_identity",
    "modeller_best_template_coverage",

    # Activity label columns should not be predictors
    "acc_activity_label_types",
    "acc_activity_source_files",
    "acc_activity_source_sheets",
    "acc_km_value_mean",
    "acc_kcat_km_value_mean",
    "acc_experimental_pi_mean",
    "acc_ph_stability_width_max",
    "acc_activity_substrates",
    "acc_activity_commentaries",
    "org_activity_label_types",
    "org_activity_source_files",
    "org_activity_source_sheets",
    "org_km_value_mean",
    "org_kcat_km_value_mean",
    "org_experimental_pi_mean",
    "org_ph_stability_width_max",
    "org_activity_substrates",
    "org_activity_commentaries",
    "activity_label_types",
    "activity_source_files",
    "activity_source_sheets",
    "activity_substrates",
    "activity_commentaries",
    "km_value_mean",
    "kcat_km_value_mean",
    "experimental_pi_mean",
    "ph_stability_width_max",
    "has_km_excel",
    "has_kcat_km_excel",
    "has_experimental_pi_excel",
    "has_ph_stability_excel",
    "has_excel_activity_label_accession",
    "has_excel_activity_label_organism",
    "has_excel_activity_label",
    "has_activity_evidence",
    "activity_match_confidence",

    # BRENDA raw labels should not be predictors in the main stability model
    "brenda_temperature_optimum",
    "brenda_temperature_range",
    "brenda_temperature_stability",
    "brenda_ph_optimum",
    "brenda_ph_range",
    "has_brenda_temp",
    "has_brenda_ph",
    "has_any_brenda",
]


def rmse(y_true, y_pred):
    return mean_squared_error(y_true, y_pred) ** 0.5


def clean_bool_like(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = (
                df[c]
                .astype(str)
                .replace({
                    "True": "true",
                    "False": "false",
                    "nan": "",
                    "NaN": "",
                    "None": "",
                })
            )
    return df


def build_preprocessor(numeric_features, categorical_features):
    numeric_transformer = Pipeline(steps=[
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
    ])

    categorical_transformer = Pipeline(steps=[
        ("imputer", SimpleImputer(strategy="most_frequent")),
        ("onehot", OneHotEncoder(handle_unknown="ignore")),
    ])

    return ColumnTransformer(
        transformers=[
            ("num", numeric_transformer, numeric_features),
            ("cat", categorical_transformer, categorical_features),
        ],
        remainder="drop",
    )


def get_feature_names(preprocessor, numeric_features, categorical_features):
    names = []

    names.extend(numeric_features)

    try:
        ohe = preprocessor.named_transformers_["cat"].named_steps["onehot"]
        cat_names = list(ohe.get_feature_names_out(categorical_features))
        names.extend(cat_names)
    except Exception:
        pass

    return names


def main():
    if not INFILE.exists():
        raise FileNotFoundError(f"Missing input file: {INFILE}")

    df = pd.read_csv(INFILE)

    if TARGET not in df.columns:
        raise ValueError(f"Target column not found: {TARGET}")

    # Keep only rows with structural stability target
    df[TARGET] = pd.to_numeric(df[TARGET], errors="coerce")
    model_df = df[df[TARGET].notna()].copy()

    # Select available, non-excluded features
    numeric_features = [
        c for c in NUMERIC_FEATURES
        if c in model_df.columns and c not in EXCLUDE_FROM_FEATURES
    ]

    categorical_features = [
        c for c in CATEGORICAL_FEATURES
        if c in model_df.columns and c not in EXCLUDE_FROM_FEATURES
    ]

    if not numeric_features and not categorical_features:
        raise ValueError("No usable ML features were found.")

    for c in numeric_features:
        model_df[c] = pd.to_numeric(model_df[c], errors="coerce")

    model_df = clean_bool_like(model_df, categorical_features)

    X = model_df[numeric_features + categorical_features].copy()
    y = model_df[TARGET].copy()

    cv = KFold(n_splits=5, shuffle=True, random_state=42)

    models = {
        "RandomForestRegressor": RandomForestRegressor(
            n_estimators=500,
            random_state=42,
            min_samples_leaf=2,
            n_jobs=-1,
        ),
        "HistGradientBoostingRegressor": HistGradientBoostingRegressor(
            random_state=42,
            max_iter=500,
            learning_rate=0.05,
            l2_regularization=0.01,
        ),
    }

    summaries = []
    prediction_tables = []
    feature_importance_tables = []
    permutation_tables = []

    for model_name, regressor in models.items():
        preprocessor = build_preprocessor(numeric_features, categorical_features)

        pipe = Pipeline(steps=[
            ("preprocessor", preprocessor),
            ("model", regressor),
        ])

        # Cross-validated predictions
        y_pred = cross_val_predict(pipe, X, y, cv=cv, n_jobs=-1)

        r2 = r2_score(y, y_pred)
        mae = mean_absolute_error(y, y_pred)
        rmse_value = rmse(y, y_pred)

        summaries.append({
            "target": TARGET,
            "model": model_name,
            "n_samples": len(model_df),
            "n_numeric_features": len(numeric_features),
            "n_categorical_features": len(categorical_features),
            "numeric_features": ",".join(numeric_features),
            "categorical_features": ",".join(categorical_features),
            "excluded_template_features": "modeller_best_template_identity,modeller_best_template_coverage",
            "cv": "5-fold KFold shuffle random_state=42",
            "r2": r2,
            "mae": mae,
            "rmse": rmse_value,
        })

        base_cols = [
            "uniprot_accession",
            "organism",
            "organism_type",
            "gh_family",
            "structure_sources",
            "foldx_source_layer",
            TARGET,
        ]

        available_base_cols = [c for c in base_cols if c in model_df.columns]

        pred_table = model_df[available_base_cols].copy()
        pred_table["model"] = model_name
        pred_table["predicted_foldx_energy_per_residue"] = y_pred
        pred_table["prediction_error"] = (
            pred_table[TARGET] - pred_table["predicted_foldx_energy_per_residue"]
        )
        pred_table["abs_error"] = pred_table["prediction_error"].abs()

        prediction_tables.append(pred_table)

        # Fit full model for feature importance/permutation importance
        pipe.fit(X, y)

        fitted_preprocessor = pipe.named_steps["preprocessor"]
        feature_names = get_feature_names(
            fitted_preprocessor,
            numeric_features,
            categorical_features,
        )

        fitted_model = pipe.named_steps["model"]

        if hasattr(fitted_model, "feature_importances_"):
            importances = fitted_model.feature_importances_
            fi = pd.DataFrame({
                "model": model_name,
                "feature": feature_names[:len(importances)],
                "importance": importances,
            }).sort_values("importance", ascending=False)

            feature_importance_tables.append(fi)

        # Permutation importance on full data
        perm = permutation_importance(
            pipe,
            X,
            y,
            n_repeats=20,
            random_state=42,
            n_jobs=-1,
            scoring="r2",
        )

        perm_df = pd.DataFrame({
            "model": model_name,
            "feature": numeric_features + categorical_features,
            "permutation_importance_mean": perm.importances_mean,
            "permutation_importance_std": perm.importances_std,
        }).sort_values("permutation_importance_mean", ascending=False)

        permutation_tables.append(perm_df)

    summary_df = pd.DataFrame(summaries)
    pred_df = pd.concat(prediction_tables, ignore_index=True)

    if feature_importance_tables:
        feature_df = pd.concat(feature_importance_tables, ignore_index=True)
    else:
        feature_df = pd.DataFrame(columns=["model", "feature", "importance"])

    perm_all = pd.concat(permutation_tables, ignore_index=True)

    summary_df.to_csv(SUMMARY_OUT, index=False)
    pred_df.to_csv(PRED_OUT, index=False)
    feature_df.to_csv(FEATURE_OUT, index=False)
    perm_all.to_csv(PERM_OUT, index=False)

    # Grouped prediction summaries
    group_cols = [
        c for c in ["model", "organism_type", "gh_family"]
        if c in pred_df.columns
    ]

    if all(c in pred_df.columns for c in ["model", "organism_type", "gh_family"]):
        group_summary = (
            pred_df.groupby(["model", "organism_type", "gh_family"])
            .agg(
                n=("uniprot_accession", "count"),
                mean_observed=(TARGET, "mean"),
                median_observed=(TARGET, "median"),
                mean_predicted=("predicted_foldx_energy_per_residue", "mean"),
                median_predicted=("predicted_foldx_energy_per_residue", "median"),
                mean_abs_error=("abs_error", "mean"),
                median_abs_error=("abs_error", "median"),
            )
            .reset_index()
        )
    else:
        group_summary = pd.DataFrame()

    group_summary.to_csv(GROUP_OUT, index=False)

    audit = []
    audit.append("Structural Thermostability ML Audit - No Template Features")
    audit.append("=" * 80)
    audit.append(f"Input file: {INFILE}")
    audit.append(f"Target: {TARGET}")
    audit.append("")
    audit.append("Template-quality features excluded:")
    audit.append("- modeller_best_template_identity")
    audit.append("- modeller_best_template_coverage")
    audit.append("")
    audit.append(f"Rows in input master: {len(df)}")
    audit.append(f"Rows used for ML: {len(model_df)}")

    if "uniprot_accession" in model_df.columns:
        audit.append(f"Unique accessions used: {model_df['uniprot_accession'].nunique()}")

    audit.append("")
    audit.append("Numeric features used:")
    for c in numeric_features:
        audit.append(f"- {c}")

    audit.append("")
    audit.append("Categorical features used:")
    for c in categorical_features:
        audit.append(f"- {c}")

    audit.append("")
    audit.append("Target summary:")
    audit.append(str(model_df[TARGET].describe()))

    audit.append("")
    audit.append("Model performance summary:")
    audit.append(summary_df.round(4).to_string(index=False))

    audit.append("")
    audit.append("Group prediction summary:")
    if len(group_summary) > 0:
        audit.append(group_summary.round(4).to_string(index=False))
    else:
        audit.append("Group summary could not be generated because grouping columns were missing.")

    audit.append("")
    audit.append("Top feature importance:")
    if len(feature_df) > 0:
        audit.append(feature_df.head(30).round(6).to_string(index=False))
    else:
        audit.append("No model-native feature importance available.")

    audit.append("")
    audit.append("Top permutation importance:")
    audit.append(perm_all.head(30).round(6).to_string(index=False))

    audit.append("")
    audit.append("Best predicted matches by absolute error:")
    audit.append(
        pred_df.sort_values("abs_error")
        .head(30)
        .round(6)
        .to_string(index=False)
    )

    audit.append("")
    audit.append("Largest prediction errors:")
    audit.append(
        pred_df.sort_values("abs_error", ascending=False)
        .head(30)
        .round(6)
        .to_string(index=False)
    )

    AUDIT_OUT.write_text("\n".join(audit))

    print(f"Saved: {SUMMARY_OUT}")
    print(f"Saved: {PRED_OUT}")
    print(f"Saved: {FEATURE_OUT}")
    print(f"Saved: {PERM_OUT}")
    print(f"Saved: {GROUP_OUT}")
    print(f"Saved: {AUDIT_OUT}")
    print("")
    print("Model performance summary:")
    print(summary_df.round(4).to_string(index=False))
    print("")
    print("Top permutation importance:")
    print(perm_all.head(20).round(6).to_string(index=False))


if __name__ == "__main__":
    main()
