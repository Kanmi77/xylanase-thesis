#!/usr/bin/env python3
"""
Regenerate all thesis visualizations with explanatory legends, readable labels,
topic-based folders, and a priority_13 folder.

Project:
In Silico Structural Thermostability Assessment of Bacterial and Fungal
Thermostable Xylanases through Bioinformatics and Machine Learning predictions.

Main output:
  results/visualization_thesis/

Subdirectories:
  01_data_curation
  02_sequence_features
  03_structural_features
  04_foldx_stability
  05_docking
  06_mutation_stability
  07_machine_learning
  08_md
  09_tmalign
  10_integrated_ranking
  priority_13
"""

from pathlib import Path
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

BASE = Path.home() / "xylanase-thesis"
OUT = BASE / "results" / "visualization_thesis"

DIRS = {
    "curation": OUT / "01_data_curation",
    "sequence": OUT / "02_sequence_features",
    "structure": OUT / "03_structural_features",
    "foldx": OUT / "04_foldx_stability",
    "docking": OUT / "05_docking",
    "mutation": OUT / "06_mutation_stability",
    "ml": OUT / "07_machine_learning",
    "md": OUT / "08_md",
    "tmalign": OUT / "09_tmalign",
    "ranking": OUT / "10_integrated_ranking",
    "priority": OUT / "priority_13",
}

for d in DIRS.values():
    d.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})


# ============================================================
# General helpers
# ============================================================

def savefig(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def read_csv(path):
    path = Path(path)
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except Exception as e:
        print(f"[WARN] Could not read {path}: {e}")
        return None


def first_existing(paths):
    for p in paths:
        p = Path(p)
        if p.exists():
            return p
    return None


def standardize_columns(df):
    if df is None:
        return None
    df = df.copy()
    df.columns = [c.strip() for c in df.columns]
    return df


def numeric(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def add_note(text):
    plt.figtext(
        0.01,
        -0.02,
        text,
        ha="left",
        fontsize=8,
        wrap=True
    )


def clean_label(x):
    x = str(x)
    x = x.replace("_", " ")
    x = re.sub(r"\s+", " ", x)
    return x


def annotate_bars(ax, fmt="{:.0f}", fontsize=8):
    for p in ax.patches:
        h = p.get_height()
        if pd.notna(h):
            ax.annotate(
                fmt.format(h),
                (p.get_x() + p.get_width() / 2, h),
                ha="center",
                va="bottom",
                fontsize=fontsize,
                xytext=(0, 3),
                textcoords="offset points"
            )


def plot_grouped_bar(df, index_col, column_col, value_col, title, xlabel, ylabel, outpath, note):
    pivot = df.pivot_table(index=index_col, columns=column_col, values=value_col, aggfunc="sum").fillna(0)
    ax = pivot.plot(kind="bar", figsize=(11, 7))
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(title=clean_label(column_col), frameon=True)
    annotate_bars(ax)
    add_note(note)
    savefig(outpath)


def plot_heatmap(matrix, title, xlabel, ylabel, outpath, cbar_label, note, fmt="{:.0f}", legend_text=None):
    plt.figure(figsize=(max(8, matrix.shape[1] * 1.1), max(5, matrix.shape[0] * 0.6)))
    im = plt.imshow(matrix.values, aspect="auto")

    plt.xticks(range(matrix.shape[1]), matrix.columns, rotation=45, ha="right")
    plt.yticks(range(matrix.shape[0]), matrix.index)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix.iloc[i, j]
            if pd.notna(val):
                try:
                    label = fmt.format(val)
                except Exception:
                    label = str(val)
                plt.text(j, i, label, ha="center", va="center", fontsize=8)

    cbar = plt.colorbar(im)
    cbar.set_label(cbar_label)

    if legend_text:
        plt.legend(handles=[Patch(label=legend_text)], loc="upper right", frameon=True)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    add_note(note)
    savefig(outpath)


def make_boxplot(df, value_col, group_col, title, xlabel, ylabel, outpath, note):
    temp = df[[value_col, group_col]].dropna().copy()
    groups = list(temp[group_col].astype(str).unique())
    data = [temp[temp[group_col].astype(str) == g][value_col].dropna() for g in groups]

    plt.figure(figsize=(10, 6))
    plt.boxplot(data, labels=groups, showmeans=True)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="", label="Mean marker shown in each box"),
            Line2D([0], [0], linestyle="-", label="Box = interquartile range; whiskers = spread")
        ],
        frameon=True,
        loc="best"
    )
    add_note(note)
    savefig(outpath)


def make_scatter_by_group(df, x_col, y_col, group_col, title, xlabel, ylabel, outpath, note):
    temp = df[[x_col, y_col, group_col]].dropna().copy()
    temp[x_col] = pd.to_numeric(temp[x_col], errors="coerce")
    temp[y_col] = pd.to_numeric(temp[y_col], errors="coerce")
    temp = temp.dropna()

    plt.figure(figsize=(10, 7))
    for group, sub in temp.groupby(group_col):
        plt.scatter(sub[x_col], sub[y_col], alpha=0.75, s=45, label=f"{group} (n={len(sub)})")

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(title=clean_label(group_col), frameon=True)
    add_note(note)
    savefig(outpath)


# ============================================================
# Load available project datasets
# ============================================================

master_path = first_existing([
    BASE / "data/curated/xylanase_thesis_master_final_v5_activity_labels.csv",
    BASE / "data/curated/xylanase_master_all_curated.csv",
    BASE / "data/curated/xylanase_master_uniprot.csv",
])

structured_path = first_existing([
    BASE / "data/curated/xylanase_structured_subset_with_foldx_norm.csv",
    BASE / "data/curated/xylanase_structured_subset_with_foldx.csv",
    BASE / "data/curated/xylanase_structured_subset.csv",
])

top15_path = first_existing([
    BASE / "results/ranking/top15_candidates_with_docking.csv",
    BASE / "report_evidence_pack/results/top15_candidates_with_docking.csv",
    BASE / "top15_candidates_with_docking.csv",
])

top15_final_path = first_existing([
    BASE / "results/ranking/top15_final_ranked.csv",
    BASE / "top15_final_ranked.csv",
])

foldx_ddg_path = first_existing([
    BASE / "results/foldx_clean/tier2_ddg_ranked_annotated.csv",
    BASE / "results/foldx_clean/tier2_ddg_ranked.csv",
])

docking_scores_path = first_existing([
    BASE / "docking_tier2_all_best_mutants_reformatted/tier2_reformatted_docking_scores.csv",
    BASE / "docking_tier2_all_best_mutants/tier2_all_best_mutant_docking_scores.csv",
    BASE / "results/docking/vina_top15_scores.csv",
])

docking_comp_path = first_existing([
    BASE / "docking_tier2_all_best_mutants_reformatted/tier2_reformatted_wt_mutant_docking_comparison.csv",
    BASE / "docking_tier2_all_best_mutants/tier2_all_best_mutant_docking_comparison.csv",
])

docking_group_path = first_existing([
    BASE / "docking_tier2_all_best_mutants_reformatted/tier2_reformatted_docking_group_summary.csv",
])

ml_predictions_path = first_existing([
    BASE / "results/ml/structural_thermostability_ml_no_template_predictions.csv",
    BASE / "results/ml/ml_cv_predictions.csv",
])

ml_feature_path = first_existing([
    BASE / "results/ml/structural_thermostability_ml_no_template_feature_importance.csv",
    BASE / "results/ml/feature_importance_all.csv",
])

ml_perm_path = first_existing([
    BASE / "results/ml/structural_thermostability_ml_no_template_permutation_importance.csv",
])

ml_group_path = first_existing([
    BASE / "results/ml/structural_thermostability_ml_no_template_group_summary.csv",
])

tm_quality_path = first_existing([
    BASE / "results/reports/tmalign_best_reference_quality_classes.csv",
])

tm_by_gh_path = first_existing([
    BASE / "results/reports/tmalign_best_reference_by_gh_family.csv",
])

tm_by_org_gh_path = first_existing([
    BASE / "results/reports/tmalign_best_reference_by_organism_gh.csv",
])

master = standardize_columns(read_csv(master_path)) if master_path else None
structured = standardize_columns(read_csv(structured_path)) if structured_path else None
top15 = standardize_columns(read_csv(top15_path)) if top15_path else None
top15_final = standardize_columns(read_csv(top15_final_path)) if top15_final_path else None
foldx_ddg = standardize_columns(read_csv(foldx_ddg_path)) if foldx_ddg_path else None
docking_scores = standardize_columns(read_csv(docking_scores_path)) if docking_scores_path else None
docking_comp = standardize_columns(read_csv(docking_comp_path)) if docking_comp_path else None
docking_group = standardize_columns(read_csv(docking_group_path)) if docking_group_path else None
ml_predictions = standardize_columns(read_csv(ml_predictions_path)) if ml_predictions_path else None
ml_feature = standardize_columns(read_csv(ml_feature_path)) if ml_feature_path else None
ml_perm = standardize_columns(read_csv(ml_perm_path)) if ml_perm_path else None
ml_group = standardize_columns(read_csv(ml_group_path)) if ml_group_path else None
tm_quality = standardize_columns(read_csv(tm_quality_path)) if tm_quality_path else None
tm_by_gh = standardize_columns(read_csv(tm_by_gh_path)) if tm_by_gh_path else None
tm_by_org_gh = standardize_columns(read_csv(tm_by_org_gh_path)) if tm_by_org_gh_path else None


# ============================================================
# 01 DATA CURATION VISUALS
# ============================================================

if master is not None:
    if {"organism_type", "gh_family"}.issubset(master.columns):
        tmp = master.copy()
        tmp["organism_type"] = tmp["organism_type"].fillna("unknown")
        tmp["gh_family"] = tmp["gh_family"].fillna("unknown")

        counts = tmp.groupby(["organism_type", "gh_family"]).size().reset_index(name="count")

        plot_grouped_bar(
            counts,
            "organism_type",
            "gh_family",
            "count",
            "Curated xylanase dataset by organism type and GH family",
            "Organism type",
            "Number of proteins",
            DIRS["curation"] / "dataset_by_organism_type_and_gh_family.png",
            "This plot checks whether bacterial/fungal and GH10/GH11 groups are represented after curation."
        )

        pivot = counts.pivot(index="organism_type", columns="gh_family", values="count").fillna(0)
        plot_heatmap(
            pivot,
            "Dataset composition heatmap",
            "GH family",
            "Organism type",
            DIRS["curation"] / "organism_type_by_gh_family_heatmap.png",
            "Number of proteins",
            "Each tile shows the number of proteins in one organism type and GH family group.",
            fmt="{:.0f}",
            legend_text="Darker tiles indicate larger groups."
        )

    if "length" in master.columns:
        master = numeric(master, ["length"])
        if "organism_type" in master.columns:
            make_boxplot(
                master,
                "length",
                "organism_type",
                "Sequence length distribution by organism type",
                "Organism type",
                "Sequence length / amino acids",
                DIRS["curation"] / "sequence_length_by_organism_type.png",
                "This figure shows whether bacterial and fungal xylanases differ in sequence length distribution."
            )

        if "gh_family" in master.columns:
            make_boxplot(
                master,
                "length",
                "gh_family",
                "Sequence length distribution by GH family",
                "GH family",
                "Sequence length / amino acids",
                DIRS["curation"] / "sequence_length_by_gh_family.png",
                "GH10 and GH11 enzymes often differ in size; this plot visualizes that structural-family difference."
            )

    if "structure_sources" in master.columns:
        ss = master["structure_sources"].fillna("unknown").value_counts().reset_index()
        ss.columns = ["structure_sources", "count"]

        plt.figure(figsize=(10, 6))
        ax = plt.bar(ss["structure_sources"], ss["count"], label="Structure source category")
        plt.title("Availability of structure sources in the curated dataset")
        plt.xlabel("Structure source")
        plt.ylabel("Number of proteins")
        plt.xticks(rotation=45, ha="right")
        plt.legend(frameon=True)
        add_note("This figure separates proteins with experimental structures from those requiring model-based structures.")
        savefig(DIRS["curation"] / "structure_source_distribution.png")


# ============================================================
# 02 SEQUENCE FEATURE VISUALS
# ============================================================

seq_source = master if master is not None else structured

if seq_source is not None:
    sequence_features = [
        "sequence_length",
        "length",
        "molecular_weight",
        "aromaticity",
        "instability_index",
        "isoelectric_point",
        "predicted_pI",
        "gravy",
        "helix_fraction",
        "turn_fraction",
        "sheet_fraction",
        "nterm_hydrophobic_count_5_25",
        "nterm_positive_count_1_10",
    ]

    available_seq = [c for c in sequence_features if c in seq_source.columns]
    seq_source = numeric(seq_source, available_seq)

    for feature in available_seq:
        if "organism_type" in seq_source.columns:
            make_boxplot(
                seq_source,
                feature,
                "organism_type",
                f"{clean_label(feature).title()} by organism type",
                "Organism type",
                clean_label(feature),
                DIRS["sequence"] / f"{feature}_by_organism_type.png",
                f"This plot compares {clean_label(feature)} between organism groups."
            )

        if "gh_family" in seq_source.columns:
            make_boxplot(
                seq_source,
                feature,
                "gh_family",
                f"{clean_label(feature).title()} by GH family",
                "GH family",
                clean_label(feature),
                DIRS["sequence"] / f"{feature}_by_gh_family.png",
                f"This plot compares {clean_label(feature)} between GH families."
            )

    corr_cols = [c for c in available_seq if seq_source[c].notna().sum() > 5]
    if len(corr_cols) >= 2:
        corr = seq_source[corr_cols].corr()
        plot_heatmap(
            corr,
            "Correlation heatmap of sequence-derived features",
            "Feature",
            "Feature",
            DIRS["sequence"] / "sequence_feature_correlation_heatmap.png",
            "Pearson correlation",
            "Positive values indicate direct relationships; negative values indicate inverse relationships.",
            fmt="{:.2f}",
            legend_text="Values are Pearson correlation coefficients."
        )


# ============================================================
# 03 STRUCTURAL FEATURE VISUALS
# ============================================================

struct_source = structured if structured is not None else master

if struct_source is not None:
    structural_features = [
        "chain_length",
        "hbond_proxy_count",
        "salt_bridge_count",
        "disulfide_count",
        "sasa_total",
        "sasa_per_residue",
        "sasa_per_res",
        "hbond_per_res",
        "salt_bridge_per_res",
        "disulfide_per_res",
        "foldx_energy_per_residue",
    ]

    available_struct = [c for c in structural_features if c in struct_source.columns]
    struct_source = numeric(struct_source, available_struct)

    for feature in available_struct:
        if "organism_type" in struct_source.columns:
            make_boxplot(
                struct_source,
                feature,
                "organism_type",
                f"{clean_label(feature).title()} by organism type",
                "Organism type",
                clean_label(feature),
                DIRS["structure"] / f"{feature}_by_organism_type.png",
                f"This structural descriptor is compared across bacterial and fungal groups."
            )

        if "gh_family" in struct_source.columns:
            make_boxplot(
                struct_source,
                feature,
                "gh_family",
                f"{clean_label(feature).title()} by GH family",
                "GH family",
                clean_label(feature),
                DIRS["structure"] / f"{feature}_by_gh_family.png",
                f"This structural descriptor is compared across GH10 and GH11 xylanases."
            )

    corr_cols = [c for c in available_struct if struct_source[c].notna().sum() > 5]
    if len(corr_cols) >= 2:
        corr = struct_source[corr_cols].corr()
        plot_heatmap(
            corr,
            "Correlation heatmap of structural stability descriptors",
            "Structural feature",
            "Structural feature",
            DIRS["structure"] / "structural_feature_correlation_heatmap.png",
            "Pearson correlation",
            "This heatmap highlights relationships between FoldX, hydrogen-bond, salt-bridge, disulfide, and SASA descriptors.",
            fmt="{:.2f}",
            legend_text="Values are Pearson correlation coefficients."
        )

    if {"hbond_proxy_count", "foldx_energy_per_residue", "organism_type"}.issubset(struct_source.columns):
        make_scatter_by_group(
            struct_source,
            "hbond_proxy_count",
            "foldx_energy_per_residue",
            "organism_type",
            "Hydrogen-bond proxy count versus FoldX energy per residue",
            "Hydrogen-bond proxy count",
            "FoldX energy per residue",
            DIRS["structure"] / "hbond_proxy_vs_foldx_energy_by_organism.png",
            "This plot tests whether hydrogen-bond-rich structures show different FoldX stability energy profiles."
        )

    if {"salt_bridge_count", "foldx_energy_per_residue", "organism_type"}.issubset(struct_source.columns):
        make_scatter_by_group(
            struct_source,
            "salt_bridge_count",
            "foldx_energy_per_residue",
            "organism_type",
            "Salt-bridge count versus FoldX energy per residue",
            "Salt-bridge count",
            "FoldX energy per residue",
            DIRS["structure"] / "salt_bridge_vs_foldx_energy_by_organism.png",
            "This plot evaluates salt bridges as structural contributors to predicted stability."
        )


# ============================================================
# 04 FOLDX STABILITY VISUALS
# ============================================================

foldx_source = struct_source if struct_source is not None else top15

if foldx_source is not None and "foldx_energy_per_residue" in foldx_source.columns:
    foldx_source = numeric(foldx_source, ["foldx_energy_per_residue", "foldx_wt_total_energy"])

    if "organism_type" in foldx_source.columns:
        make_boxplot(
            foldx_source,
            "foldx_energy_per_residue",
            "organism_type",
            "FoldX energy per residue by organism type",
            "Organism type",
            "FoldX energy per residue",
            DIRS["foldx"] / "foldx_energy_per_residue_by_organism_type.png",
            "Lower FoldX energy per residue is interpreted as a more favorable structural stability proxy."
        )

    if "gh_family" in foldx_source.columns:
        make_boxplot(
            foldx_source,
            "foldx_energy_per_residue",
            "gh_family",
            "FoldX energy per residue by GH family",
            "GH family",
            "FoldX energy per residue",
            DIRS["foldx"] / "foldx_energy_per_residue_by_gh_family.png",
            "This figure compares structural stability proxy values between GH10 and GH11 xylanases."
        )

    if {"organism_type", "gh_family"}.issubset(foldx_source.columns):
        grouped = foldx_source.groupby(["organism_type", "gh_family"])["foldx_energy_per_residue"].mean().reset_index()
        plot_grouped_bar(
            grouped,
            "organism_type",
            "gh_family",
            "foldx_energy_per_residue",
            "Mean FoldX energy per residue by organism type and GH family",
            "Organism type",
            "Mean FoldX energy per residue",
            DIRS["foldx"] / "mean_foldx_energy_per_residue_by_organism_gh.png",
            "This grouped bar chart summarizes the stability proxy across the main comparison groups."
        )

if foldx_ddg is not None:
    ddg_col = "ddg" if "ddg" in foldx_ddg.columns else "foldx_ddg" if "foldx_ddg" in foldx_ddg.columns else None

    if ddg_col:
        foldx_ddg = numeric(foldx_ddg, [ddg_col])

        plt.figure(figsize=(10, 6))
        plt.hist(foldx_ddg[ddg_col].dropna(), bins=35, alpha=0.75, label=f"All mutations (n={foldx_ddg[ddg_col].notna().sum()})")
        plt.axvline(0, linestyle="--", linewidth=1, label="ddG = 0 threshold")
        plt.title("Distribution of FoldX-predicted mutation effects")
        plt.xlabel("FoldX ddG")
        plt.ylabel("Frequency")
        plt.legend(frameon=True)
        add_note("Negative ddG values indicate predicted stabilizing mutations; positive values indicate destabilizing mutations.")
        savefig(DIRS["foldx"] / "foldx_ddg_distribution_all_mutations.png")

        if "gh_family" in foldx_ddg.columns:
            plt.figure(figsize=(10, 6))
            for fam, sub in foldx_ddg.groupby("gh_family"):
                plt.hist(sub[ddg_col].dropna(), bins=30, alpha=0.55, label=f"{fam} (n={sub[ddg_col].notna().sum()})")
            plt.axvline(0, linestyle="--", linewidth=1, label="ddG = 0 threshold")
            plt.title("FoldX ddG distribution by GH family")
            plt.xlabel("FoldX ddG")
            plt.ylabel("Frequency")
            plt.legend(title="GH family", frameon=True)
            add_note("This plot compares predicted stabilizing mutation distributions between GH families.")
            savefig(DIRS["foldx"] / "foldx_ddg_distribution_by_gh_family.png")


# ============================================================
# 05 DOCKING VISUALS
# ============================================================

if docking_scores is not None:
    docking_scores = numeric(docking_scores, ["best_binding_energy", "mean_top3_binding_energy", "vina_best_xylobiose", "vina_best_xylotriose"])

    if {"ligand", "state", "best_binding_energy"}.issubset(docking_scores.columns):
        plt.figure(figsize=(10, 6))
        for (lig, state), sub in docking_scores.groupby(["ligand", "state"]):
            plt.hist(sub["best_binding_energy"].dropna(), bins=25, alpha=0.55, label=f"{lig} / {state} (n={len(sub)})")
        plt.title("Best docking binding energy distribution by ligand and protein state")
        plt.xlabel("Best binding energy")
        plt.ylabel("Frequency")
        plt.legend(title="Ligand / state", frameon=True)
        add_note("More negative binding energy indicates stronger predicted ligand binding.")
        savefig(DIRS["docking"] / "best_binding_energy_by_ligand_and_state.png")

    if {"ligand", "best_binding_energy"}.issubset(docking_scores.columns):
        make_boxplot(
            docking_scores,
            "best_binding_energy",
            "ligand",
            "Best docking binding energy by ligand",
            "Ligand",
            "Best binding energy",
            DIRS["docking"] / "best_binding_energy_by_ligand.png",
            "This plot compares xylobiose and xylotriose docking strength across the docked systems."
        )

    if {"state", "best_binding_energy"}.issubset(docking_scores.columns):
        make_boxplot(
            docking_scores,
            "best_binding_energy",
            "state",
            "Best docking binding energy by WT/mutant state",
            "Protein state",
            "Best binding energy",
            DIRS["docking"] / "best_binding_energy_by_state.png",
            "This plot compares ligand binding distributions between wild-type and mutant structures."
        )

if docking_comp is not None:
    docking_comp = numeric(
        docking_comp,
        [
            "wt_binding_energy",
            "mut_binding_energy",
            "delta_binding_mut_minus_wt",
            "wt_mean_top3_binding",
            "mut_mean_top3_binding",
            "delta_top3_mut_minus_wt",
            "foldx_ddg",
            "ddg",
            "foldx_energy_per_residue",
            "sasa_per_res",
            "hbond_per_res",
            "disulfide_per_res",
        ]
    )

    if {"wt_binding_energy", "mut_binding_energy"}.issubset(docking_comp.columns):
        plt.figure(figsize=(10, 6))
        plt.hist(docking_comp["wt_binding_energy"].dropna(), bins=30, alpha=0.60, label=f"WT (n={docking_comp['wt_binding_energy'].notna().sum()})")
        plt.hist(docking_comp["mut_binding_energy"].dropna(), bins=30, alpha=0.60, label=f"Mutant (n={docking_comp['mut_binding_energy'].notna().sum()})")
        plt.title("WT and mutant docking binding energy distributions")
        plt.xlabel("Binding energy")
        plt.ylabel("Frequency")
        plt.legend(title="Protein state", frameon=True)
        add_note("A shift toward more negative values indicates stronger predicted ligand binding.")
        savefig(DIRS["docking"] / "wt_mutant_binding_energy_distribution.png")

    if {"delta_binding_mut_minus_wt", "functional_integrity"}.issubset(docking_comp.columns):
        make_boxplot(
            docking_comp,
            "delta_binding_mut_minus_wt",
            "functional_integrity",
            "Mutation effect on docking binding energy by functional integrity class",
            "Functional integrity class",
            "Δ binding energy, mutant - WT",
            DIRS["docking"] / "delta_binding_by_functional_integrity.png",
            "Negative Δ binding means the mutant binds more strongly than WT; positive values indicate weaker mutant binding."
        )

    if {"foldx_ddg", "delta_binding_mut_minus_wt", "functional_integrity"}.issubset(docking_comp.columns):
        make_scatter_by_group(
            docking_comp,
            "foldx_ddg",
            "delta_binding_mut_minus_wt",
            "functional_integrity",
            "FoldX ddG versus docking binding change",
            "FoldX ddG",
            "Δ binding energy, mutant - WT",
            DIRS["docking"] / "foldx_ddg_vs_delta_binding_by_functional_integrity.png",
            "This figure evaluates whether stabilizing mutations also retain or improve ligand binding."
        )

    if {"ligand", "delta_binding_mut_minus_wt"}.issubset(docking_comp.columns):
        make_boxplot(
            docking_comp,
            "delta_binding_mut_minus_wt",
            "ligand",
            "Mutation effect on binding energy by ligand",
            "Ligand",
            "Δ binding energy, mutant - WT",
            DIRS["docking"] / "delta_binding_by_ligand.png",
            "This plot compares mutation effects for xylobiose and xylotriose docking."
        )

if docking_group is not None:
    docking_group = numeric(
        docking_group,
        [
            "mean_wt_binding",
            "mean_mut_binding",
            "mean_delta_binding",
            "retained_or_improved",
            "improved",
            "weakened",
            "mean_foldx_ddg",
            "retained_or_improved_fraction",
            "proteins",
            "comparisons",
        ]
    )

    if {"organism_type", "gh_family", "ligand", "retained_or_improved_fraction"}.issubset(docking_group.columns):
        df = docking_group.copy()
        df["group"] = df["organism_type"].astype(str) + " / " + df["gh_family"].astype(str) + " / " + df["ligand"].astype(str)

        plt.figure(figsize=(13, 6))
        plt.bar(df["group"], df["retained_or_improved_fraction"], label="Retained or improved fraction")
        plt.xticks(rotation=45, ha="right")
        plt.ylim(0, 1.1)
        plt.title("Functional retention fraction after mutation by organism type, GH family, and ligand")
        plt.xlabel("Group")
        plt.ylabel("Retained or improved fraction")
        plt.legend(frameon=True)
        for i, v in enumerate(df["retained_or_improved_fraction"]):
            if pd.notna(v):
                plt.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)
        add_note("Values close to 1 indicate that most mutations retained or improved docking-based functional binding.")
        savefig(DIRS["docking"] / "functional_retention_fraction_by_group.png")


# ============================================================
# 06 MUTATION STABILITY VISUALS
# ============================================================

if foldx_ddg is not None:
    ddg_col = "ddg" if "ddg" in foldx_ddg.columns else "foldx_ddg" if "foldx_ddg" in foldx_ddg.columns else None

    if ddg_col:
        if "organism_type" in foldx_ddg.columns:
            make_boxplot(
                foldx_ddg,
                ddg_col,
                "organism_type",
                "FoldX mutation ddG by organism type",
                "Organism type",
                "FoldX ddG",
                DIRS["mutation"] / "foldx_ddg_by_organism_type.png",
                "Negative values indicate predicted stabilizing mutations."
            )

        if "gh_family" in foldx_ddg.columns:
            make_boxplot(
                foldx_ddg,
                ddg_col,
                "gh_family",
                "FoldX mutation ddG by GH family",
                "GH family",
                "FoldX ddG",
                DIRS["mutation"] / "foldx_ddg_by_gh_family.png",
                "This figure compares predicted mutational stabilization across GH families."
            )

        if {"organism_type", "gh_family"}.issubset(foldx_ddg.columns):
            summary = foldx_ddg.copy()
            summary["stabilizing"] = summary[ddg_col] < 0
            grouped = summary.groupby(["organism_type", "gh_family"])["stabilizing"].mean().reset_index()
            grouped["stabilizing_fraction"] = grouped["stabilizing"]

            plot_grouped_bar(
                grouped,
                "organism_type",
                "gh_family",
                "stabilizing_fraction",
                "Fraction of predicted stabilizing mutations by organism type and GH family",
                "Organism type",
                "Fraction of mutations with ddG < 0",
                DIRS["mutation"] / "stabilizing_mutation_fraction_by_organism_gh.png",
                "This plot shows the proportion of tested mutations predicted to stabilize each group."
            )


# ============================================================
# 07 MACHINE LEARNING VISUALS
# ============================================================

if ml_predictions is not None:
    possible_obs = [
        "foldx_energy_per_residue",
        "observed_foldx_energy_per_residue",
        "y_true",
        "observed",
        "actual",
    ]
    possible_pred = [
        "predicted_foldx_energy_per_residue",
        "y_pred",
        "predicted",
        "prediction",
    ]

    obs_col = next((c for c in possible_obs if c in ml_predictions.columns), None)
    pred_col = next((c for c in possible_pred if c in ml_predictions.columns), None)

    if obs_col is None:
        for c in ml_predictions.columns:
            if "observed" in c.lower() or "actual" in c.lower() or "true" in c.lower():
                obs_col = c
                break

    if pred_col is None:
        for c in ml_predictions.columns:
            if "pred" in c.lower():
                pred_col = c
                break

    if obs_col and pred_col:
        ml_predictions = numeric(ml_predictions, [obs_col, pred_col])
        df = ml_predictions[[obs_col, pred_col]].dropna()

        plt.figure(figsize=(8, 7))
        plt.scatter(df[obs_col], df[pred_col], alpha=0.70, s=45, label=f"Predicted samples (n={len(df)})")
        minv = min(df[obs_col].min(), df[pred_col].min())
        maxv = max(df[obs_col].max(), df[pred_col].max())
        plt.plot([minv, maxv], [minv, maxv], linestyle="--", label="Ideal prediction line")

        corr = df[[obs_col, pred_col]].corr().iloc[0, 1]
        plt.title("Observed versus predicted FoldX energy per residue")
        plt.xlabel("Observed FoldX energy per residue")
        plt.ylabel("Predicted FoldX energy per residue")
        plt.legend(frameon=True)
        add_note(f"Pearson correlation between observed and predicted values = {corr:.3f}.")
        savefig(DIRS["ml"] / "observed_vs_predicted_foldx_energy_per_residue.png")

        df["residual"] = df[pred_col] - df[obs_col]
        plt.figure(figsize=(10, 6))
        plt.hist(df["residual"], bins=30, alpha=0.75, label=f"Residuals (n={len(df)})")
        plt.axvline(0, linestyle="--", linewidth=1, label="Zero error")
        plt.title("Machine-learning prediction residual distribution")
        plt.xlabel("Prediction residual, predicted - observed")
        plt.ylabel("Frequency")
        plt.legend(frameon=True)
        add_note("Residuals close to zero indicate better prediction accuracy.")
        savefig(DIRS["ml"] / "ml_prediction_residual_distribution.png")

if ml_feature is not None:
    feature_col = next((c for c in ml_feature.columns if "feature" in c.lower()), None)
    importance_col = next((c for c in ml_feature.columns if "importance" in c.lower()), None)

    if feature_col and importance_col:
        ml_feature = numeric(ml_feature, [importance_col])
        df = ml_feature[[feature_col, importance_col]].dropna().sort_values(importance_col, ascending=False).head(20)

        plt.figure(figsize=(11, 8))
        plt.barh(df[feature_col][::-1], df[importance_col][::-1], label="Feature importance")
        plt.title("Top machine-learning feature importance values")
        plt.xlabel("Importance")
        plt.ylabel("Feature")
        plt.legend(frameon=True)
        add_note("Higher values indicate stronger contribution to the ML prediction of the structural thermostability proxy.")
        savefig(DIRS["ml"] / "ml_top_feature_importance.png")

if ml_perm is not None:
    feature_col = next((c for c in ml_perm.columns if "feature" in c.lower()), None)
    perm_col = next((c for c in ml_perm.columns if "permutation_importance_mean" in c.lower()), None)

    if feature_col and perm_col:
        ml_perm = numeric(ml_perm, [perm_col])
        df = ml_perm[[feature_col, perm_col]].dropna().sort_values(perm_col, ascending=False).head(20)

        plt.figure(figsize=(11, 8))
        plt.barh(df[feature_col][::-1], df[perm_col][::-1], label="Permutation importance")
        plt.title("Top permutation importance features")
        plt.xlabel("Mean permutation importance")
        plt.ylabel("Feature")
        plt.legend(frameon=True)
        add_note("Permutation importance estimates how much model performance decreases when each feature is disrupted.")
        savefig(DIRS["ml"] / "ml_top_permutation_importance.png")

if ml_group is not None and {"model", "organism_type", "gh_family", "mean_abs_error"}.issubset(ml_group.columns):
    ml_group = numeric(ml_group, ["mean_abs_error", "median_abs_error", "n"])
    df = ml_group.copy()
    df["group"] = df["organism_type"].astype(str) + " / " + df["gh_family"].astype(str)

    plt.figure(figsize=(12, 6))
    for model, sub in df.groupby("model"):
        plt.plot(sub["group"], sub["mean_abs_error"], marker="o", label=model)

    plt.xticks(rotation=45, ha="right")
    plt.title("ML mean absolute error by organism type and GH family")
    plt.xlabel("Group")
    plt.ylabel("Mean absolute error")
    plt.legend(title="Model", frameon=True)
    add_note("Lower mean absolute error indicates better group-level prediction performance.")
    savefig(DIRS["ml"] / "ml_group_prediction_error_by_model.png")


# ============================================================
# 08 MD STATUS VISUALS
# ============================================================

def classify_md(row):
    text = str(row.get("status_text", "")).lower()
    if "done" in text or "complete" in text or "completed" in text:
        return "Completed"
    if row.get("has_md_gro", False) and row.get("has_md_xtc", False):
        return "Trajectory output present"
    if row.get("has_md_log", False) or row.get("has_md_cpt", False) or row.get("has_md_xtc", False):
        return "Started / partial"
    return "Not started"


md_roots = [
    BASE / "md_tier2_wt_mutant_compact/systems",
    BASE / "md_final8_apo/systems",
    BASE / "md_tier2_mini/systems",
    BASE / "md/systems",
]

md_rows = []

for root in md_roots:
    if not root.exists():
        continue

    for system_dir in sorted(root.iterdir()):
        if not system_dir.is_dir():
            continue

        for temp_dir in sorted(system_dir.iterdir()):
            if not temp_dir.is_dir():
                continue

            if not re.search(r"\d+K$", temp_dir.name):
                continue

            status_file = temp_dir / "status.txt"
            row = {
                "md_root": root.name,
                "system": system_dir.name,
                "temperature": temp_dir.name,
                "has_md_log": (temp_dir / "md.log").exists(),
                "has_md_xtc": (temp_dir / "md.xtc").exists(),
                "has_md_gro": (temp_dir / "md.gro").exists(),
                "has_md_cpt": (temp_dir / "md.cpt").exists(),
                "status_text": status_file.read_text(errors="ignore").strip() if status_file.exists() else "",
            }
            row["status"] = classify_md(row)
            md_rows.append(row)

if md_rows:
    md = pd.DataFrame(md_rows)

    status_order = {
        "Not started": 0,
        "Started / partial": 1,
        "Trajectory output present": 2,
        "Completed": 3,
    }

    md["status_num"] = md["status"].map(status_order)

    pivot = md.pivot_table(index="system", columns="temperature", values="status_num", aggfunc="max").fillna(0)
    label_pivot = md.pivot_table(index="system", columns="temperature", values="status", aggfunc="first")

    plt.figure(figsize=(11, max(5, 0.45 * len(pivot))))
    im = plt.imshow(pivot.values, aspect="auto", vmin=0, vmax=3)

    plt.xticks(range(len(pivot.columns)), pivot.columns)
    plt.yticks(range(len(pivot.index)), pivot.index, fontsize=7)

    short = {
        "Not started": "NS",
        "Started / partial": "SP",
        "Trajectory output present": "TO",
        "Completed": "C",
    }

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = label_pivot.iloc[i, j]
            if pd.notna(val):
                plt.text(j, i, short.get(val, ""), ha="center", va="center", fontsize=7)

    cbar = plt.colorbar(im)
    cbar.set_label("MD completion status")
    cbar.set_ticks([0, 1, 2, 3])
    cbar.set_ticklabels(["Not started", "Partial", "Trajectory output", "Completed"])

    plt.legend(
        handles=[
            Patch(label="NS = not started"),
            Patch(label="SP = started or partial files present"),
            Patch(label="TO = trajectory output present"),
            Patch(label="C = completed"),
        ],
        title="Status legend",
        frameon=True,
        loc="upper right"
    )

    plt.title("Molecular dynamics simulation status overview")
    plt.xlabel("Temperature")
    plt.ylabel("MD system")
    add_note("This plot summarizes current MD execution status across detected systems and simulation temperatures.")
    savefig(DIRS["md"] / "md_simulation_status_heatmap.png")

    status_counts = md.groupby(["temperature", "status"]).size().reset_index(name="count")
    plot_grouped_bar(
        status_counts,
        "temperature",
        "status",
        "count",
        "MD status count by temperature",
        "Temperature",
        "Number of systems",
        DIRS["md"] / "md_status_count_by_temperature.png",
        "This grouped plot shows how many simulations are completed, partial, or still missing at each temperature."
    )


# ============================================================
# 09 TM-ALIGN VISUALS
# ============================================================

if tm_quality is not None:
    rename_map = {}
    for c in tm_quality.columns:
        lc = c.lower()
        if lc in ["tm_score", "best_tm_score"]:
            rename_map[c] = "tm_score_best"
        if lc in ["rmsd_best", "best_rmsd"]:
            rename_map[c] = "rmsd"
    tm_quality = tm_quality.rename(columns=rename_map)

    if {"tm_score_best", "rmsd"}.issubset(tm_quality.columns):
        tm_quality = numeric(tm_quality, ["tm_score_best", "rmsd"])

        class_col = None
        for c in tm_quality.columns:
            if "quality" in c.lower() or "class" in c.lower():
                class_col = c
                break

        if class_col is None:
            class_col = "quality_class"
            tm_quality[class_col] = "all_models"

        plt.figure(figsize=(10, 7))
        for cls, sub in tm_quality.dropna(subset=["tm_score_best", "rmsd"]).groupby(class_col):
            plt.scatter(sub["rmsd"], sub["tm_score_best"], alpha=0.75, s=45, label=f"{cls} (n={len(sub)})")

        plt.axhline(0.95, linestyle="--", linewidth=1, label="TM-score 0.95 reference")
        plt.axvline(1.0, linestyle="--", linewidth=1, label="RMSD 1.0 Å reference")
        plt.title("TM-align structural validation: TM-score versus RMSD")
        plt.xlabel("RMSD / Å")
        plt.ylabel("Best TM-score")
        plt.legend(title="Model quality", frameon=True)
        add_note("High TM-score and low RMSD indicate stronger agreement between generated models and structural references.")
        savefig(DIRS["tmalign"] / "tm_score_vs_rmsd_with_quality_legend.png")

    if class_col and class_col in tm_quality.columns:
        counts = tm_quality[class_col].fillna("unknown").value_counts().reset_index()
        counts.columns = [class_col, "count"]

        plt.figure(figsize=(10, 6))
        plt.bar(counts[class_col], counts["count"], label="Number of structures")
        plt.title("TM-align model quality class distribution")
        plt.xlabel("Quality class")
        plt.ylabel("Number of models")
        plt.xticks(rotation=45, ha="right")
        plt.legend(frameon=True)
        for i, v in enumerate(counts["count"]):
            plt.text(i, v + 0.5, str(v), ha="center", fontsize=8)
        add_note("This plot summarizes how many structures fall into each TM-align quality category.")
        savefig(DIRS["tmalign"] / "tmalign_quality_class_distribution.png")

if tm_by_gh is not None:
    possible_tm = next((c for c in tm_by_gh.columns if "tm" in c.lower() and "score" in c.lower()), None)
    if possible_tm and "gh_family" in tm_by_gh.columns:
        tm_by_gh = numeric(tm_by_gh, [possible_tm])
        make_boxplot(
            tm_by_gh,
            possible_tm,
            "gh_family",
            "TM-align score distribution by GH family",
            "GH family",
            "TM-score",
            DIRS["tmalign"] / "tmalign_score_by_gh_family.png",
            "This plot checks whether model-reference agreement differs between GH10 and GH11 proteins."
        )


# ============================================================
# 10 INTEGRATED RANKING VISUALS
# ============================================================

ranking_source = top15_final if top15_final is not None else top15

if ranking_source is not None:
    score_cols = [
        "final_score",
        "composite_score",
        "composite_rank",
        "foldx_energy_per_residue",
        "foldx_wt_total_energy",
        "salt_bridge_count",
        "hbond_proxy_count",
        "hbond_per_res",
        "salt_bridge_per_res",
        "disulfide_count",
        "sasa_total",
        "sasa_per_residue",
        "sasa_per_res",
        "vina_best_xylobiose",
        "vina_best_xylotriose",
        "vina_best_min",
        "vina_best_mean",
        "vina_xylobiose_score",
        "vina_xylotriose_score",
    ]

    available_scores = [c for c in score_cols if c in ranking_source.columns]
    ranking_source = numeric(ranking_source, available_scores)

    id_col = None
    for c in ["pdb_tag", "uniprot_accession", "protein", "pdb_id"]:
        if c in ranking_source.columns:
            id_col = c
            break

    sort_col = None
    for c in ["final_score", "composite_rank", "composite_score", "foldx_energy_per_residue"]:
        if c in ranking_source.columns:
            sort_col = c
            break

    if id_col and sort_col:
        df = ranking_source.copy()
        ascending = True if "rank" in sort_col or "energy" in sort_col else False
        df = df.sort_values(sort_col, ascending=ascending).head(20)

        plt.figure(figsize=(12, 7))
        plt.barh(df[id_col].astype(str)[::-1], df[sort_col][::-1], label=clean_label(sort_col))
        plt.title(f"Top ranked candidate xylanases by {clean_label(sort_col)}")
        plt.xlabel(clean_label(sort_col))
        plt.ylabel("Candidate")
        plt.legend(frameon=True)
        add_note("This plot summarizes the leading candidates according to the selected integrated ranking metric.")
        savefig(DIRS["ranking"] / f"top_candidates_by_{sort_col}.png")

    if id_col and len(available_scores) >= 3:
        df = ranking_source.copy()
        if sort_col:
            ascending = True if "rank" in sort_col or "energy" in sort_col else False
            df = df.sort_values(sort_col, ascending=ascending).head(20)
        else:
            df = df.head(20)

        heat = df[[id_col] + available_scores].set_index(id_col)
        heat = heat.apply(pd.to_numeric, errors="coerce")

        # Min-max normalize column-wise
        norm = heat.copy()
        for c in norm.columns:
            col = norm[c]
            if col.notna().sum() > 1 and col.max() != col.min():
                norm[c] = (col - col.min()) / (col.max() - col.min())
            else:
                norm[c] = 0.0

        plot_heatmap(
            norm,
            "Integrated candidate prioritization heatmap",
            "Normalized feature",
            "Candidate",
            DIRS["ranking"] / "integrated_candidate_prioritization_heatmap.png",
            "Min-max normalized value",
            "Each feature is normalized independently to compare candidates across different scoring scales.",
            fmt="{:.2f}",
            legend_text="Values are normalized within each feature column."
        )

    if {"organism_type", "gh_family"}.issubset(ranking_source.columns) and sort_col:
        grouped = ranking_source.groupby(["organism_type", "gh_family"])[sort_col].mean().reset_index()
        plot_grouped_bar(
            grouped,
            "organism_type",
            "gh_family",
            sort_col,
            f"Mean {clean_label(sort_col)} by organism type and GH family",
            "Organism type",
            f"Mean {clean_label(sort_col)}",
            DIRS["ranking"] / f"mean_{sort_col}_by_organism_gh.png",
            "This plot compares candidate prioritization patterns across the main biological groups."
        )


# ============================================================
# Copy selected plots into priority_13
# ============================================================

priority_map = [
    (DIRS["curation"] / "dataset_by_organism_type_and_gh_family.png", "01_dataset_curation_overview.png"),
    (DIRS["curation"] / "organism_type_by_gh_family_heatmap.png", "02_organism_type_by_gh_family_heatmap.png"),
    (DIRS["tmalign"] / "tm_score_vs_rmsd_with_quality_legend.png", "03_tm_score_vs_rmsd.png"),
    (DIRS["foldx"] / "mean_foldx_energy_per_residue_by_organism_gh.png", "04_foldx_energy_per_residue_by_organism_gh.png"),
    (DIRS["structure"] / "structural_feature_correlation_heatmap.png", "05_structural_feature_correlation_heatmap.png"),
    (DIRS["foldx"] / "foldx_ddg_distribution_all_mutations.png", "06_foldx_ddg_distribution.png"),
    (DIRS["mutation"] / "stabilizing_mutation_fraction_by_organism_gh.png", "07_mutation_effect_by_organism_gh.png"),
    (DIRS["docking"] / "wt_mutant_binding_energy_distribution.png", "08_wt_mutant_binding_energy_distribution.png"),
    (DIRS["docking"] / "foldx_ddg_vs_delta_binding_by_functional_integrity.png", "09_foldx_ddg_vs_delta_binding.png"),
    (DIRS["ml"] / "observed_vs_predicted_foldx_energy_per_residue.png", "10_observed_vs_predicted_foldx_energy.png"),
    (DIRS["ml"] / "ml_top_feature_importance.png", "11_ml_feature_importance.png"),
    (DIRS["md"] / "md_simulation_status_heatmap.png", "12_md_simulation_status_heatmap.png"),
    (DIRS["ranking"] / "integrated_candidate_prioritization_heatmap.png", "13_integrated_candidate_prioritization_heatmap.png"),
]

for src, name in priority_map:
    if src.exists():
        dst = DIRS["priority"] / name
        dst.write_bytes(src.read_bytes())
        print(f"Priority copy: {dst}")
    else:
        print(f"[WARN] Priority source missing: {src}")


# ============================================================
# Audit report
# ============================================================

audit = OUT / "visualization_regeneration_audit.txt"

with open(audit, "w") as f:
    f.write("Thesis visualization regeneration audit\n")
    f.write("=" * 80 + "\n\n")

    f.write("Input files detected:\n")
    for label, path in [
        ("master", master_path),
        ("structured", structured_path),
        ("top15", top15_path),
        ("top15_final", top15_final_path),
        ("foldx_ddg", foldx_ddg_path),
        ("docking_scores", docking_scores_path),
        ("docking_comparison", docking_comp_path),
        ("docking_group", docking_group_path),
        ("ml_predictions", ml_predictions_path),
        ("ml_feature", ml_feature_path),
        ("ml_permutation", ml_perm_path),
        ("ml_group", ml_group_path),
        ("tmalign_quality", tm_quality_path),
        ("tmalign_by_gh", tm_by_gh_path),
        ("tmalign_by_organism_gh", tm_by_org_gh_path),
    ]:
        f.write(f"- {label}: {path if path else 'MISSING'}\n")

    f.write("\nGenerated PNG files:\n")
    for p in sorted(OUT.rglob("*.png")):
        f.write(f"- {p.relative_to(OUT)}\n")

print(f"Saved audit: {audit}")
print("\nAll thesis visualizations regenerated with legends and explanatory annotations.")
print(f"Main output directory: {OUT}")
print(f"Priority 13 directory: {DIRS['priority']}")
