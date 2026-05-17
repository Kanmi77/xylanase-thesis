#!/usr/bin/env python3

from pathlib import Path
import shutil
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

BASE = Path.home() / "xylanase-thesis"

VIS = BASE / "results/visualization_thesis"
PRIORITY = VIS / "priority_13"

DIRS = {
    "data_curation": VIS / "01_data_acquisition_curation",
    "sequence": VIS / "02_sequence_features",
    "phylogeny": VIS / "03_phylogeny_conservation",
    "tmalign": VIS / "04_structural_modelling_tmalign",
    "structural": VIS / "05_structural_features_foldx",
    "mutation": VIS / "06_mutation_screening",
    "docking": VIS / "07_docking",
    "ml": VIS / "08_machine_learning",
    "md": VIS / "09_molecular_dynamics",
    "integration": VIS / "10_integrated_ranking",
}

for d in DIRS.values():
    d.mkdir(parents=True, exist_ok=True)
PRIORITY.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------

def savefig(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def read_csv_safe(path):
    path = Path(path)
    if not path.exists():
        print(f"[MISSING] {path}")
        return None
    try:
        return pd.read_csv(path)
    except Exception as e:
        print(f"[ERROR READING] {path}: {e}")
        return None


def numeric(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def bar_plot(series, title, ylabel, outpath, rotation=0):
    series = series.dropna()
    plt.figure(figsize=(8, 5))
    plt.bar(series.index.astype(str), series.values)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=rotation, ha="right" if rotation else "center")
    for i, v in enumerate(series.values):
        plt.text(i, v, str(int(v)), ha="center", va="bottom", fontsize=9)
    savefig(outpath)


def grouped_bar(df, x, y, group, title, ylabel, outpath, rotation=30):
    pivot = df.pivot(index=x, columns=group, values=y).fillna(0)
    ax = pivot.plot(kind="bar", figsize=(10, 6))
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    plt.xticks(rotation=rotation, ha="right")
    plt.legend(title=group, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(outpath)


def heatmap_from_table(table, title, outpath, fmt="{:.0f}"):
    fig, ax = plt.subplots(figsize=(8, 5))
    data = table.values.astype(float)
    im = ax.imshow(data, aspect="auto")
    ax.set_xticks(np.arange(table.shape[1]))
    ax.set_yticks(np.arange(table.shape[0]))
    ax.set_xticklabels(table.columns.astype(str), rotation=30, ha="right")
    ax.set_yticklabels(table.index.astype(str))
    ax.set_title(title)

    for i in range(table.shape[0]):
        for j in range(table.shape[1]):
            val = data[i, j]
            text = fmt.format(val) if not np.isnan(val) else ""
            ax.text(j, i, text, ha="center", va="center", fontsize=9)

    fig.colorbar(im, ax=ax)
    savefig(outpath)


def boxplot_by_group(df, value, group, title, ylabel, outpath, rotation=30):
    if value not in df.columns or group not in df.columns:
        return

    data = []
    labels = []

    for name, sub in df.groupby(group):
        vals = pd.to_numeric(sub[value], errors="coerce").dropna()
        if len(vals) > 0:
            data.append(vals.values)
            labels.append(str(name))

    if not data:
        return

    plt.figure(figsize=(10, 6))
    plt.boxplot(data, labels=labels, showfliers=False)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=rotation, ha="right")
    savefig(outpath)


def boxplot_by_two_groups(df, value, group1, group2, title, ylabel, outpath):
    if value not in df.columns or group1 not in df.columns or group2 not in df.columns:
        return

    df = df.copy()
    df["combined_group"] = df[group1].astype(str) + " " + df[group2].astype(str)

    data = []
    labels = []

    for name, sub in df.groupby("combined_group"):
        vals = pd.to_numeric(sub[value], errors="coerce").dropna()
        if len(vals) > 0:
            data.append(vals.values)
            labels.append(name)

    if not data:
        return

    plt.figure(figsize=(12, 6))
    plt.boxplot(data, labels=labels, showfliers=False)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=35, ha="right")
    savefig(outpath)


def histogram(df, value, title, xlabel, outpath, vline=None):
    if value not in df.columns:
        return

    vals = pd.to_numeric(df[value], errors="coerce").dropna()
    if len(vals) == 0:
        return

    plt.figure(figsize=(8, 5))
    plt.hist(vals, bins=40)
    if vline is not None:
        plt.axvline(vline, linestyle="--")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    savefig(outpath)


def scatter(df, x, y, title, xlabel, ylabel, outpath, color_col=None, diagonal=False):
    if x not in df.columns or y not in df.columns:
        return

    plot_df = df.copy()
    plot_df[x] = pd.to_numeric(plot_df[x], errors="coerce")
    plot_df[y] = pd.to_numeric(plot_df[y], errors="coerce")
    plot_df = plot_df.dropna(subset=[x, y])

    if len(plot_df) == 0:
        return

    plt.figure(figsize=(8, 6))

    if color_col and color_col in plot_df.columns:
        groups = plot_df[color_col].fillna("unknown").astype(str)
        for g in sorted(groups.unique()):
            sub = plot_df[groups == g]
            plt.scatter(sub[x], sub[y], s=25, alpha=0.75, label=g)
        plt.legend(title=color_col, bbox_to_anchor=(1.02, 1), loc="upper left")
    else:
        plt.scatter(plot_df[x], plot_df[y], s=25, alpha=0.75)

    if diagonal:
        lo = min(plot_df[x].min(), plot_df[y].min())
        hi = max(plot_df[x].max(), plot_df[y].max())
        plt.plot([lo, hi], [lo, hi], linestyle="--")
        plt.xlim(lo, hi)
        plt.ylim(lo, hi)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    savefig(outpath)


def correlation_heatmap(df, cols, title, outpath):
    cols = [c for c in cols if c in df.columns]
    if len(cols) < 2:
        return

    tmp = df[cols].copy()
    for c in cols:
        tmp[c] = pd.to_numeric(tmp[c], errors="coerce")

    corr = tmp.corr()
    if corr.empty:
        return

    fig, ax = plt.subplots(figsize=(11, 9))
    im = ax.imshow(corr.values, vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(corr.columns)))
    ax.set_yticks(np.arange(len(corr.index)))
    ax.set_xticklabels(corr.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(corr.index, fontsize=8)
    ax.set_title(title)

    for i in range(corr.shape[0]):
        for j in range(corr.shape[1]):
            val = corr.iloc[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)

    fig.colorbar(im, ax=ax)
    savefig(outpath)


def horizontal_bar(df, label_col, value_col, title, xlabel, outpath, top_n=20, ascending=True):
    if label_col not in df.columns or value_col not in df.columns:
        return

    tmp = df[[label_col, value_col]].copy()
    tmp[value_col] = pd.to_numeric(tmp[value_col], errors="coerce")
    tmp = tmp.dropna(subset=[value_col])
    tmp = tmp.sort_values(value_col, ascending=ascending).head(top_n)

    if len(tmp) == 0:
        return

    plt.figure(figsize=(10, max(5, 0.35 * len(tmp))))
    y = np.arange(len(tmp))
    plt.barh(y, tmp[value_col].values)
    plt.yticks(y, tmp[label_col].astype(str).values)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel(xlabel)
    savefig(outpath)


def copy_priority(src_path, priority_name):
    src_path = Path(src_path)
    if src_path.exists():
        dst = PRIORITY / priority_name
        shutil.copy2(src_path, dst)
        print(f"Priority copy: {dst}")


# ---------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------

MASTER = BASE / "data/curated/xylanase_thesis_master_final_v5_activity_labels.csv"
TMALIGN = BASE / "results/reports/tmalign_best_reference_per_model.csv"
TMALIGN_ALT = BASE / "results/reports/tmalign_best_reference_summary.csv"

MUT_ANNOTATED = BASE / "results/foldx_clean/tier2_ddg_ranked_annotated.csv"
MUT_BY_GROUP = BASE / "results/foldx_clean/tier2_ddg_summary_by_group.csv"

DOCK_SCORES = BASE / "docking_tier2_all_best_mutants_reformatted/tier2_reformatted_docking_scores.csv"
DOCK_COMP = BASE / "docking_tier2_all_best_mutants_reformatted/tier2_reformatted_wt_mutant_docking_comparison.csv"
DOCK_GROUP = BASE / "docking_tier2_all_best_mutants_reformatted/tier2_reformatted_docking_group_summary.csv"
DOCK_RETAINED = BASE / "docking_tier2_all_best_mutants_reformatted/tier2_reformatted_top_functionally_retained_mutants.csv"

ML_SUMMARY = BASE / "results/ml/structural_thermostability_ml_no_template_summary.csv"
ML_PRED = BASE / "results/ml/structural_thermostability_ml_no_template_predictions.csv"
ML_FEATURE = BASE / "results/ml/structural_thermostability_ml_no_template_feature_importance.csv"
ML_PERM = BASE / "results/ml/structural_thermostability_ml_no_template_permutation_importance.csv"
ML_GROUP = BASE / "results/ml/structural_thermostability_ml_no_template_group_summary.csv"

master = read_csv_safe(MASTER)
tmalign = read_csv_safe(TMALIGN)
if tmalign is None:
    tmalign = read_csv_safe(TMALIGN_ALT)

mut = read_csv_safe(MUT_ANNOTATED)
mut_group = read_csv_safe(MUT_BY_GROUP)

dock_scores = read_csv_safe(DOCK_SCORES)
dock_comp = read_csv_safe(DOCK_COMP)
dock_group = read_csv_safe(DOCK_GROUP)
dock_retained = read_csv_safe(DOCK_RETAINED)

ml_summary = read_csv_safe(ML_SUMMARY)
ml_pred = read_csv_safe(ML_PRED)
ml_feature = read_csv_safe(ML_FEATURE)
ml_perm = read_csv_safe(ML_PERM)
ml_group = read_csv_safe(ML_GROUP)


# ---------------------------------------------------------------------
# 01 Data acquisition and curation
# ---------------------------------------------------------------------

if master is not None:
    master = numeric(master, [
        "sequence_length", "molecular_weight", "aromaticity", "instability_index",
        "isoelectric_point", "gravy", "helix_fraction", "turn_fraction",
        "sheet_fraction", "chain_length", "foldx_wt_total_energy",
        "foldx_energy_per_residue", "hbond_proxy_count", "salt_bridge_count",
        "disulfide_count", "sasa_total", "hbond_per_res", "salt_bridge_per_res",
        "disulfide_per_res", "sasa_per_res"
    ])

    # 1. Dataset construction flowchart as bar representation
    stage_counts = pd.Series({
        "Initial UniProt\nxylanase entries": 36337,
        "Curated GH10/GH11\nworking dataset": len(master),
        "Modelled\nstructures": int((master.get("structure_sources", pd.Series()).fillna("") == "modeller").sum()) if "structure_sources" in master else 631,
        "PDB\nstructures": int((master.get("structure_sources", pd.Series()).fillna("") == "pdb").sum()) if "structure_sources" in master else 48,
        "No structure": int(master["structure_sources"].isna().sum()) if "structure_sources" in master else 13,
        "Activity\nevidence": int(master["has_activity_evidence"].astype(str).str.lower().eq("true").sum()) if "has_activity_evidence" in master else 187,
    })
    bar_plot(
        stage_counts,
        "Overview of xylanase data acquisition and curation",
        "Number of entries",
        DIRS["data_curation"] / "figure_4_01_dataset_curation_overview.png",
        rotation=20
    )

    # 2. Organism distribution
    if "organism_type" in master.columns:
        bar_plot(
            master["organism_type"].fillna("unknown").value_counts(),
            "Distribution of curated GH10/GH11 xylanases by organism type",
            "Number of proteins",
            DIRS["data_curation"] / "figure_4_02_organism_type_distribution.png"
        )

    # 3. GH distribution
    if "gh_family" in master.columns:
        bar_plot(
            master["gh_family"].fillna("unknown").value_counts(),
            "Distribution of curated xylanases by GH family",
            "Number of proteins",
            DIRS["data_curation"] / "figure_4_03_gh_family_distribution.png"
        )

    # 4. organism x GH heatmap
    if {"organism_type", "gh_family"}.issubset(master.columns):
        table = pd.crosstab(master["organism_type"], master["gh_family"])
        heatmap_from_table(
            table,
            "Cross-distribution of curated xylanases by organism type and GH family",
            DIRS["data_curation"] / "figure_4_04_organism_type_by_gh_family_heatmap.png"
        )

    # 5. Activity evidence coverage
    if "activity_match_confidence" in master.columns:
        bar_plot(
            master["activity_match_confidence"].fillna("none").value_counts(),
            "Activity-label evidence coverage in the curated xylanase dataset",
            "Number of proteins",
            DIRS["data_curation"] / "figure_4_05_activity_evidence_confidence.png",
            rotation=30
        )

    # 6. Activity evidence by organism x GH
    if {"organism_type", "gh_family", "has_activity_evidence"}.issubset(master.columns):
        tmp = master.copy()
        tmp["has_activity_evidence"] = tmp["has_activity_evidence"].astype(str).str.lower().map({"true": "evidence", "false": "no evidence"}).fillna("no evidence")
        tmp["group"] = tmp["organism_type"].astype(str) + " " + tmp["gh_family"].astype(str)
        pivot = pd.crosstab(tmp["group"], tmp["has_activity_evidence"])
        ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 6))
        ax.set_title("Distribution of activity evidence across organism type and GH family")
        ax.set_ylabel("Number of proteins")
        ax.set_xlabel("")
        plt.xticks(rotation=35, ha="right")
        plt.legend(title="Activity evidence", bbox_to_anchor=(1.02, 1), loc="upper left")
        savefig(DIRS["data_curation"] / "figure_4_06_activity_evidence_by_organism_gh.png")

    # Structure source distribution
    if "structure_sources" in master.columns:
        bar_plot(
            master["structure_sources"].fillna("no_structure").value_counts(),
            "Structure-source distribution in the curated xylanase dataset",
            "Number of proteins",
            DIRS["data_curation"] / "figure_4_14_structure_source_distribution.png"
        )


# ---------------------------------------------------------------------
# 02 Sequence features
# ---------------------------------------------------------------------

if master is not None:
    boxplot_by_group(
        master, "sequence_length", "gh_family",
        "Sequence length distribution of GH10 and GH11 xylanases",
        "Sequence length (aa)",
        DIRS["sequence"] / "figure_4_07_sequence_length_by_gh_family.png"
    )

    boxplot_by_two_groups(
        master, "sequence_length", "organism_type", "gh_family",
        "Sequence length variation across organism type and GH family",
        "Sequence length (aa)",
        DIRS["sequence"] / "figure_4_08_sequence_length_by_organism_gh.png"
    )

    seq_features = [
        "molecular_weight", "aromaticity", "instability_index",
        "isoelectric_point", "gravy", "helix_fraction",
        "turn_fraction", "sheet_fraction"
    ]

    for feat in seq_features:
        boxplot_by_two_groups(
            master, feat, "organism_type", "gh_family",
            f"{feat.replace('_', ' ').title()} across curated xylanase groups",
            feat.replace("_", " ").title(),
            DIRS["sequence"] / f"figure_4_09_{feat}_by_organism_gh.png"
        )

    correlation_heatmap(
        master,
        ["sequence_length"] + seq_features,
        "Correlation matrix of sequence-derived physicochemical features",
        DIRS["sequence"] / "figure_4_10_sequence_feature_correlation_heatmap.png"
    )


# ---------------------------------------------------------------------
# 03 Phylogeny/conservation placeholder folder
# ---------------------------------------------------------------------

phylo_note = DIRS["phylogeny"] / "README_phylogeny_visuals.txt"
phylo_note.write_text(
    "Place exported phylogenetic tree images and sequence-logo figures here.\n\n"
    "Recommended figures:\n"
    "1. Bacterial GH10 tree\n"
    "2. Bacterial GH11 tree\n"
    "3. Fungal GH10 tree\n"
    "4. Fungal GH11 tree\n"
    "5. Combined tree annotated with top FoldX/mutation/docking/MD candidates\n"
    "6. GH10 and GH11 sequence-logo conservation plots if generated from alignments\n"
)
print(f"Saved: {phylo_note}")


# ---------------------------------------------------------------------
# 04 Structural modelling and TM-align
# ---------------------------------------------------------------------

if tmalign is not None:
    tmalign = numeric(tmalign, ["tm_score_best", "rmsd", "seq_id_aligned"])

    histogram(
        tmalign, "tm_score_best",
        "Distribution of best-reference TM-align scores",
        "TM-score",
        DIRS["tmalign"] / "figure_4_15_tmalign_tm_score_distribution.png"
    )

    histogram(
        tmalign, "rmsd",
        "RMSD distribution from TM-align best-reference comparisons",
        "RMSD (Å)",
        DIRS["tmalign"] / "figure_4_16_tmalign_rmsd_distribution.png"
    )

    scatter(
        tmalign, "tm_score_best", "rmsd",
        "Relationship between TM-score and RMSD in model-reference structural alignment",
        "TM-score",
        "RMSD (Å)",
        DIRS["tmalign"] / "figure_4_17_tm_score_vs_rmsd.png",
        color_col="model_gh_family" if "model_gh_family" in tmalign.columns else None
    )

    if "tmalign_quality_class" in tmalign.columns:
        bar_plot(
            tmalign["tmalign_quality_class"].fillna("unknown").value_counts(),
            "TM-align structural quality classification of modelled xylanases",
            "Number of models",
            DIRS["tmalign"] / "figure_4_18_tmalign_quality_class_distribution.png"
        )

    if {"model_organism_type", "model_gh_family", "tmalign_quality_class"}.issubset(tmalign.columns):
        tmp = tmalign.copy()
        tmp["group"] = tmp["model_organism_type"].astype(str) + " " + tmp["model_gh_family"].astype(str)
        pivot = pd.crosstab(tmp["group"], tmp["tmalign_quality_class"])
        ax = pivot.plot(kind="bar", stacked=True, figsize=(11, 6))
        ax.set_title("TM-align quality classes across organism type and GH family")
        ax.set_ylabel("Number of models")
        ax.set_xlabel("")
        plt.xticks(rotation=35, ha="right")
        plt.legend(title="Quality class", bbox_to_anchor=(1.02, 1), loc="upper left")
        savefig(DIRS["tmalign"] / "figure_4_19_tmalign_quality_by_organism_gh.png")


# ---------------------------------------------------------------------
# 05 Structural features and FoldX
# ---------------------------------------------------------------------

if master is not None:
    histogram(
        master, "foldx_energy_per_residue",
        "Distribution of FoldX-normalized energy per residue",
        "FoldX energy per residue",
        DIRS["structural"] / "figure_4_20_foldx_energy_per_residue_distribution.png"
    )

    boxplot_by_two_groups(
        master, "foldx_energy_per_residue", "organism_type", "gh_family",
        "FoldX energy per residue across organism type and GH family",
        "FoldX energy per residue",
        DIRS["structural"] / "figure_4_21_foldx_energy_per_residue_by_organism_gh.png"
    )

    boxplot_by_two_groups(
        master, "hbond_proxy_count", "organism_type", "gh_family",
        "Hydrogen-bond proxy count distribution across organism type and GH family",
        "H-bond proxy count",
        DIRS["structural"] / "figure_4_22_hbond_proxy_by_organism_gh.png"
    )

    boxplot_by_two_groups(
        master, "salt_bridge_count", "organism_type", "gh_family",
        "Salt-bridge count distribution across organism type and GH family",
        "Salt bridge count",
        DIRS["structural"] / "figure_4_23_salt_bridge_count_by_organism_gh.png"
    )

    boxplot_by_two_groups(
        master, "sasa_per_res", "organism_type", "gh_family",
        "SASA per residue across organism type and GH family",
        "SASA per residue",
        DIRS["structural"] / "figure_4_24_sasa_per_residue_by_organism_gh.png"
    )

    structural_cols = [
        "foldx_energy_per_residue", "chain_length", "hbond_proxy_count",
        "salt_bridge_count", "disulfide_count", "sasa_total",
        "hbond_per_res", "salt_bridge_per_res", "disulfide_per_res",
        "sasa_per_res"
    ]

    correlation_heatmap(
        master, structural_cols,
        "Correlation matrix of structural features and FoldX-normalized energy",
        DIRS["structural"] / "figure_4_25_structural_feature_correlation_heatmap.png"
    )

    for x in ["hbond_per_res", "salt_bridge_per_res", "sasa_per_res", "chain_length"]:
        scatter(
            master, x, "foldx_energy_per_residue",
            f"Relationship between FoldX-normalized energy and {x.replace('_', ' ')}",
            x.replace("_", " ").title(),
            "FoldX energy per residue",
            DIRS["structural"] / f"figure_4_26_foldx_energy_vs_{x}.png",
            color_col="gh_family"
        )


# ---------------------------------------------------------------------
# 06 Mutation screening
# ---------------------------------------------------------------------

if mut is not None:
    mut = numeric(mut, ["ddg", "foldx_energy_per_residue", "sasa_per_res", "hbond_per_res", "disulfide_per_res"])

    mutation_counts = pd.Series({
        "Proteins tested": mut["protein"].nunique() if "protein" in mut.columns else 150,
        "Mutations evaluated": len(mut),
        "Stabilizing\nΔΔG < 0": int((mut["ddg"] < 0).sum()) if "ddg" in mut.columns else 875,
        "Destabilizing\nΔΔG > 0": int((mut["ddg"] > 0).sum()) if "ddg" in mut.columns else 1375,
    })

    bar_plot(
        mutation_counts,
        "Tier-2 mutation screening workflow and FoldX ΔΔG classification",
        "Count",
        DIRS["mutation"] / "figure_4_27_mutation_screening_workflow_counts.png",
        rotation=20
    )

    histogram(
        mut, "ddg",
        "Distribution of FoldX ΔΔG values for Tier-2 mutations",
        "FoldX ΔΔG",
        DIRS["mutation"] / "figure_4_28_foldx_ddg_distribution.png",
        vline=0
    )

    if "organism_type" in mut.columns:
        tmp = mut.copy()
        tmp["effect_class"] = np.where(tmp["ddg"] < 0, "stabilizing", "destabilizing")
        pivot = pd.crosstab(tmp["organism_type"], tmp["effect_class"])
        ax = pivot.plot(kind="bar", stacked=True, figsize=(9, 6))
        ax.set_title("Distribution of stabilizing and destabilizing mutations by organism type")
        ax.set_ylabel("Number of mutations")
        ax.set_xlabel("")
        plt.xticks(rotation=25, ha="right")
        plt.legend(title="Mutation effect", bbox_to_anchor=(1.02, 1), loc="upper left")
        savefig(DIRS["mutation"] / "figure_4_29_mutation_effect_by_organism_type.png")

    if "gh_family" in mut.columns:
        tmp = mut.copy()
        tmp["effect_class"] = np.where(tmp["ddg"] < 0, "stabilizing", "destabilizing")
        pivot = pd.crosstab(tmp["gh_family"], tmp["effect_class"])
        ax = pivot.plot(kind="bar", stacked=True, figsize=(8, 6))
        ax.set_title("Stabilizing mutation fraction by GH family")
        ax.set_ylabel("Number of mutations")
        ax.set_xlabel("")
        plt.xticks(rotation=0)
        plt.legend(title="Mutation effect", bbox_to_anchor=(1.02, 1), loc="upper left")
        savefig(DIRS["mutation"] / "figure_4_30_mutation_effect_by_gh_family.png")

    if {"organism_type", "gh_family"}.issubset(mut.columns):
        tmp = mut.copy()
        tmp["effect_class"] = np.where(tmp["ddg"] < 0, "stabilizing", "destabilizing")
        tmp["group"] = tmp["organism_type"].astype(str) + " " + tmp["gh_family"].astype(str)
        pivot = pd.crosstab(tmp["group"], tmp["effect_class"])
        ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 6))
        ax.set_title("Mutation effect distribution across organism type and GH family")
        ax.set_ylabel("Number of mutations")
        ax.set_xlabel("")
        plt.xticks(rotation=35, ha="right")
        plt.legend(title="Mutation effect", bbox_to_anchor=(1.02, 1), loc="upper left")
        savefig(DIRS["mutation"] / "figure_4_31_mutation_effect_by_organism_gh.png")

    if {"uniprot_accession", "mutation", "ddg"}.issubset(mut.columns):
        tmp = mut.copy()
        tmp["label"] = tmp["uniprot_accession"].astype(str) + " " + tmp["mutation"].astype(str)
        horizontal_bar(
            tmp, "label", "ddg",
            "Top-ranked stabilizing mutations based on FoldX ΔΔG",
            "FoldX ΔΔG",
            DIRS["mutation"] / "figure_4_32_top_stabilizing_mutations.png",
            top_n=20,
            ascending=True
        )


# ---------------------------------------------------------------------
# 07 Docking
# ---------------------------------------------------------------------

if dock_scores is not None:
    dock_scores = numeric(dock_scores, ["best_binding_energy", "mean_top3_binding_energy", "foldx_ddg", "n_modes"])

    docking_counts = pd.Series({
        "Proteins": dock_scores["protein"].nunique() if "protein" in dock_scores.columns else 150,
        "Docking records": len(dock_scores),
        "WT records": int((dock_scores["state"] == "WT").sum()) if "state" in dock_scores.columns else 300,
        "Mutant records": int((dock_scores["state"] == "MUT").sum()) if "state" in dock_scores.columns else 300,
        "Ligands": dock_scores["ligand"].nunique() if "ligand" in dock_scores.columns else 2,
        "Vina result poses": int(dock_scores["n_modes"].sum()) if "n_modes" in dock_scores.columns else 6000,
    })

    bar_plot(
        docking_counts,
        "Coverage of Tier-2 WT/mutant docking against xylobiose and xylotriose",
        "Count",
        DIRS["docking"] / "figure_4_33_tier2_docking_coverage.png",
        rotation=25
    )

    if {"state", "ligand", "best_binding_energy"}.issubset(dock_scores.columns):
        tmp = dock_scores.copy()
        tmp["state_ligand"] = tmp["state"].astype(str) + " " + tmp["ligand"].astype(str)
        boxplot_by_group(
            tmp, "best_binding_energy", "state_ligand",
            "Distribution of WT and mutant docking energies for xylobiose and xylotriose",
            "Best binding energy",
            DIRS["docking"] / "figure_4_34_wt_mutant_binding_energy_distribution.png"
        )

if dock_comp is not None:
    dock_comp = numeric(dock_comp, [
        "foldx_ddg", "wt_binding_energy", "mut_binding_energy",
        "delta_binding_mut_minus_wt", "wt_mean_top3_binding",
        "mut_mean_top3_binding", "delta_top3_mut_minus_wt", "ddg"
    ])

    histogram(
        dock_comp, "delta_binding_mut_minus_wt",
        "Distribution of docking-energy changes between mutant and wild-type receptors",
        "Δ binding energy: mutant − WT",
        DIRS["docking"] / "figure_4_35_delta_binding_distribution.png",
        vline=0
    )

    if {"foldx_ddg", "delta_binding_mut_minus_wt", "ligand"}.issubset(dock_comp.columns):
        scatter(
            dock_comp, "foldx_ddg", "delta_binding_mut_minus_wt",
            "Relationship between FoldX-predicted stabilisation and docking-energy retention",
            "FoldX ΔΔG",
            "Δ binding energy: mutant − WT",
            DIRS["docking"] / "figure_4_39_foldx_ddg_vs_delta_binding.png",
            color_col="ligand"
        )

    if {"protein", "mutation", "ligand", "delta_binding_mut_minus_wt"}.issubset(dock_comp.columns):
        tmp = dock_comp.copy()
        tmp["label"] = tmp["protein"].astype(str) + " " + tmp["mutation"].astype(str) + " " + tmp["ligand"].astype(str)

        horizontal_bar(
            tmp, "label", "delta_binding_mut_minus_wt",
            "Mutants with strongest improved docking energies relative to wild type",
            "Δ binding energy: mutant − WT",
            DIRS["docking"] / "figure_4_37_top_improved_docking_mutants.png",
            top_n=20,
            ascending=True
        )

        horizontal_bar(
            tmp, "label", "delta_binding_mut_minus_wt",
            "Mutants with weakened docking energies relative to wild type",
            "Δ binding energy: mutant − WT",
            DIRS["docking"] / "figure_4_38_top_weakened_docking_mutants.png",
            top_n=20,
            ascending=False
        )

if dock_group is not None:
    dock_group = numeric(dock_group, [
        "proteins", "comparisons", "mean_wt_binding", "mean_mut_binding",
        "mean_delta_binding", "retained_or_improved", "improved", "weakened",
        "mean_foldx_ddg", "retained_or_improved_fraction"
    ])

    if {"organism_type", "gh_family", "ligand", "retained_or_improved_fraction"}.issubset(dock_group.columns):
        tmp = dock_group.copy()
        tmp["group"] = tmp["organism_type"].astype(str) + " " + tmp["gh_family"].astype(str)
        grouped_bar(
            tmp,
            x="group",
            y="retained_or_improved_fraction",
            group="ligand",
            title="Fraction of mutants retaining or improving docking performance across groups",
            ylabel="Retained/improved fraction",
            outpath=DIRS["docking"] / "figure_4_36_retained_improved_fraction_by_group.png"
        )


# ---------------------------------------------------------------------
# 08 Machine learning
# ---------------------------------------------------------------------

if ml_summary is not None:
    ml_summary = numeric(ml_summary, ["r2", "mae", "rmse"])

    if {"model", "r2", "mae", "rmse"}.issubset(ml_summary.columns):
        for metric in ["r2", "mae", "rmse"]:
            tmp = ml_summary[["model", metric]].dropna()
            plt.figure(figsize=(8, 5))
            plt.bar(tmp["model"].astype(str), tmp[metric])
            plt.title(f"Machine-learning model comparison: {metric.upper()}")
            plt.ylabel(metric.upper())
            plt.xticks(rotation=20, ha="right")
            savefig(DIRS["ml"] / f"figure_4_40_ml_model_comparison_{metric}.png")

if ml_pred is not None:
    ml_pred = numeric(ml_pred, [
        "foldx_energy_per_residue",
        "predicted_foldx_energy_per_residue",
        "prediction_error",
        "abs_error"
    ])

    scatter(
        ml_pred,
        "foldx_energy_per_residue",
        "predicted_foldx_energy_per_residue",
        "Observed versus predicted FoldX-normalized energy",
        "Observed FoldX energy per residue",
        "Predicted FoldX energy per residue",
        DIRS["ml"] / "figure_4_41_observed_vs_predicted_foldx_energy.png",
        color_col="organism_type" if "organism_type" in ml_pred.columns else None,
        diagonal=True
    )

    histogram(
        ml_pred,
        "prediction_error",
        "Distribution of machine-learning prediction errors",
        "Prediction error",
        DIRS["ml"] / "figure_4_42_prediction_error_distribution.png",
        vline=0
    )

    if "abs_error" in ml_pred.columns:
        histogram(
            ml_pred,
            "abs_error",
            "Distribution of absolute machine-learning prediction errors",
            "Absolute error",
            DIRS["ml"] / "figure_4_42b_absolute_error_distribution.png"
        )

if ml_feature is not None:
    ml_feature = numeric(ml_feature, ["importance"])
    if {"feature", "importance"}.issubset(ml_feature.columns):
        horizontal_bar(
            ml_feature,
            "feature",
            "importance",
            "Feature importance for FoldX energy-per-residue prediction",
            "Importance",
            DIRS["ml"] / "figure_4_43_ml_feature_importance.png",
            top_n=25,
            ascending=False
        )

if ml_perm is not None:
    ml_perm = numeric(ml_perm, ["permutation_importance_mean", "permutation_importance_std"])
    if {"feature", "permutation_importance_mean"}.issubset(ml_perm.columns):
        tmp = ml_perm.sort_values("permutation_importance_mean", ascending=False).head(25)
        plt.figure(figsize=(10, 8))
        y = np.arange(len(tmp))
        xerr = tmp["permutation_importance_std"] if "permutation_importance_std" in tmp.columns else None
        plt.barh(y, tmp["permutation_importance_mean"], xerr=xerr)
        plt.yticks(y, tmp["feature"].astype(str))
        plt.gca().invert_yaxis()
        plt.title("Permutation importance of sequence and structural features")
        plt.xlabel("Permutation importance")
        savefig(DIRS["ml"] / "figure_4_44_ml_permutation_importance.png")


# ---------------------------------------------------------------------
# 09 Molecular dynamics
# ---------------------------------------------------------------------

md_base = BASE / "md_tier2_wt_mutant_compact" / "systems"
md_status_rows = []

if md_base.exists():
    for status in sorted(md_base.glob("*/*/status.txt")):
        system = status.parent.parent.name
        temp = status.parent.name
        try:
            content = status.read_text(errors="ignore").strip()
        except Exception:
            content = ""
        md_status_rows.append({
            "system": system,
            "temperature": temp,
            "status_text": content,
            "has_md_log": (status.parent / "md.log").exists(),
            "has_xtc": (status.parent / "md.xtc").exists(),
            "has_gro": (status.parent / "md.gro").exists(),
            "has_cpt": (status.parent / "md.cpt").exists(),
        })

if md_status_rows:
    md_status = pd.DataFrame(md_status_rows)
    md_status.to_csv(DIRS["md"] / "md_simulation_status_summary.csv", index=False)

    status_matrix = md_status.pivot_table(
        index="system",
        columns="temperature",
        values="has_xtc",
        aggfunc="max",
        fill_value=False
    ).astype(int)

    heatmap_from_table(
        status_matrix,
        "Current MD simulation status for selected WT and mutant systems",
        DIRS["md"] / "figure_4_45_md_simulation_status_heatmap.png"
    )
else:
    md_note = DIRS["md"] / "README_md_visuals.txt"
    md_note.write_text(
        "MD result files were not detected or MD is still running.\n\n"
        "Recommended figures after trajectory analysis:\n"
        "1. RMSD over time: WT vs mutant at 333 K and 373 K\n"
        "2. RMSF per residue: WT vs mutant at 333 K and 373 K\n"
        "3. Radius of gyration over time\n"
        "4. Intramolecular hydrogen bonds over time\n"
        "5. SASA over time\n"
    )
    print(f"Saved: {md_note}")


# ---------------------------------------------------------------------
# 10 Integrated ranking
# ---------------------------------------------------------------------

if dock_retained is not None:
    dock_retained = numeric(dock_retained, [
        "foldx_ddg", "wt_binding_energy", "mut_binding_energy",
        "delta_binding_mut_minus_wt", "wt_mean_top3_binding",
        "mut_mean_top3_binding", "delta_top3_mut_minus_wt",
        "ddg", "foldx_energy_per_residue", "sasa_per_res",
        "hbond_per_res", "disulfide_per_res"
    ])

    # Integrated candidate table
    keep_cols = [
        "protein", "organism", "organism_type", "gh_family", "mutation",
        "ligand", "foldx_ddg", "wt_binding_energy", "mut_binding_energy",
        "delta_binding_mut_minus_wt", "functional_integrity",
        "foldx_energy_per_residue", "sasa_per_res", "hbond_per_res",
        "disulfide_per_res"
    ]
    keep_cols = [c for c in keep_cols if c in dock_retained.columns]
    candidate_table = dock_retained[keep_cols].copy()
    candidate_table.to_csv(DIRS["integration"] / "table_final_integrated_candidate_summary.csv", index=False)
    print(f"Saved: {DIRS['integration'] / 'table_final_integrated_candidate_summary.csv'}")

    # Integrated heatmap using numeric columns only
    numeric_cols = [
        "foldx_ddg", "foldx_energy_per_residue", "delta_binding_mut_minus_wt",
        "sasa_per_res", "hbond_per_res", "disulfide_per_res"
    ]
    numeric_cols = [c for c in numeric_cols if c in dock_retained.columns]

    if {"protein", "mutation"}.issubset(dock_retained.columns) and numeric_cols:
        top = dock_retained.copy()
        top["candidate"] = top["protein"].astype(str) + " " + top["mutation"].astype(str) + " " + top["ligand"].astype(str)
        top = top.sort_values(["foldx_ddg", "delta_binding_mut_minus_wt"], ascending=[True, True]).head(30)
        mat = top.set_index("candidate")[numeric_cols]
        # z-score each column for heatmap comparability
        mat_z = mat.copy()
        for c in mat_z.columns:
            s = pd.to_numeric(mat_z[c], errors="coerce")
            if s.std() and not np.isnan(s.std()):
                mat_z[c] = (s - s.mean()) / s.std()
            else:
                mat_z[c] = 0

        fig, ax = plt.subplots(figsize=(10, 10))
        im = ax.imshow(mat_z.values, aspect="auto")
        ax.set_xticks(np.arange(mat_z.shape[1]))
        ax.set_yticks(np.arange(mat_z.shape[0]))
        ax.set_xticklabels(mat_z.columns, rotation=35, ha="right")
        ax.set_yticklabels(mat_z.index, fontsize=7)
        ax.set_title("Integrated prioritization matrix of selected thermostability-enhancing xylanase mutants")
        fig.colorbar(im, ax=ax)
        savefig(DIRS["integration"] / "figure_4_51_integrated_candidate_prioritization_heatmap.png")


# ---------------------------------------------------------------------
# Copy 13 priority figures
# ---------------------------------------------------------------------

priority_map = {
    DIRS["data_curation"] / "figure_4_01_dataset_curation_overview.png":
        "01_dataset_curation_overview.png",

    DIRS["data_curation"] / "figure_4_04_organism_type_by_gh_family_heatmap.png":
        "02_organism_type_by_gh_family_heatmap.png",

    DIRS["tmalign"] / "figure_4_17_tm_score_vs_rmsd.png":
        "03_tm_score_vs_rmsd.png",

    DIRS["structural"] / "figure_4_21_foldx_energy_per_residue_by_organism_gh.png":
        "04_foldx_energy_per_residue_by_organism_gh.png",

    DIRS["structural"] / "figure_4_25_structural_feature_correlation_heatmap.png":
        "05_structural_feature_correlation_heatmap.png",

    DIRS["mutation"] / "figure_4_28_foldx_ddg_distribution.png":
        "06_foldx_ddg_distribution.png",

    DIRS["mutation"] / "figure_4_31_mutation_effect_by_organism_gh.png":
        "07_mutation_effect_by_organism_gh.png",

    DIRS["docking"] / "figure_4_34_wt_mutant_binding_energy_distribution.png":
        "08_wt_mutant_binding_energy_distribution.png",

    DIRS["docking"] / "figure_4_39_foldx_ddg_vs_delta_binding.png":
        "09_foldx_ddg_vs_delta_binding.png",

    DIRS["ml"] / "figure_4_41_observed_vs_predicted_foldx_energy.png":
        "10_observed_vs_predicted_foldx_energy.png",

    DIRS["ml"] / "figure_4_43_ml_feature_importance.png":
        "11_ml_feature_importance.png",

    DIRS["md"] / "figure_4_45_md_simulation_status_heatmap.png":
        "12_md_simulation_status_heatmap.png",

    DIRS["integration"] / "figure_4_51_integrated_candidate_prioritization_heatmap.png":
        "13_integrated_candidate_prioritization_heatmap.png",
}

for src, dst_name in priority_map.items():
    copy_priority(src, dst_name)


# ---------------------------------------------------------------------
# Write visual inventory
# ---------------------------------------------------------------------

inventory = []
for png in sorted(VIS.rglob("*.png")):
    inventory.append(str(png.relative_to(VIS)))

inventory_path = VIS / "visualization_inventory.txt"
inventory_path.write_text("\n".join(inventory))
print(f"Saved: {inventory_path}")

print("\nAll visualization generation completed.")
print(f"Main visualization directory: {VIS}")
print(f"Priority 13 directory: {PRIORITY}")
