#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path.home() / "xylanase-thesis"

MASTER_IN = BASE / "data/curated/xylanase_thesis_master_final.csv"
PDB_FOLDX = BASE / "results/foldx/foldx_with_structure.csv"

MASTER_OUT = BASE / "data/curated/xylanase_thesis_master_final_v2.csv"
AUDIT_OUT = BASE / "results/reports/xylanase_thesis_master_final_v2_audit.txt"

df = pd.read_csv(MASTER_IN)
pdb = pd.read_csv(PDB_FOLDX)

df["uniprot_accession"] = df["uniprot_accession"].astype(str).str.strip()
pdb["uniprot_accession"] = pdb["uniprot_accession"].astype(str).str.strip()

# Convert numeric columns in PDB FoldX table
for col in ["foldx_wt_total_energy", "chain_length", "hbond_count", "salt_bridges", "disulfides", "sasa"]:
    if col in pdb.columns:
        pdb[col] = pd.to_numeric(pdb[col], errors="coerce")

# Use sequence_length from accession master as fallback for missing PDB chain_length
length_map = pd.to_numeric(df.set_index("uniprot_accession")["sequence_length"], errors="coerce").to_dict()
pdb["sequence_length_fallback"] = pdb["uniprot_accession"].map(length_map)

pdb["effective_chain_length"] = pdb["chain_length"]
pdb.loc[pdb["effective_chain_length"].isna(), "effective_chain_length"] = pdb.loc[
    pdb["effective_chain_length"].isna(), "sequence_length_fallback"
]

# Calculate normalized fields using effective_chain_length
pdb["foldx_energy_per_residue_calc"] = pdb["foldx_wt_total_energy"] / pdb["effective_chain_length"]
pdb["hbond_per_res_calc"] = pdb["hbond_count"] / pdb["effective_chain_length"]
pdb["disulfide_per_res_calc"] = pdb["disulfides"] / pdb["effective_chain_length"]
pdb["sasa_per_res_calc"] = pdb["sasa"] / pdb["effective_chain_length"]

# Select one representative PDB structure per accession:
# lowest FoldX energy per residue
pdb["sort_energy"] = pdb["foldx_energy_per_residue_calc"]
pdb.loc[pdb["sort_energy"].isna(), "sort_energy"] = pdb["foldx_wt_total_energy"]

pdb_best = (
    pdb.sort_values(["uniprot_accession", "sort_energy"], ascending=[True, True])
       .groupby("uniprot_accession", as_index=False)
       .head(1)
       .set_index("uniprot_accession")
)

# Ensure target columns exist
targets = [
    "foldx_structure_source",
    "foldx_structure_id",
    "foldx_input_pdb",
    "foldx_ready",
    "chain_length",
    "foldx_wt_total_energy",
    "foldx_energy_per_residue",
    "hbond_proxy_count",
    "salt_bridge_count",
    "disulfide_count",
    "sasa_total",
    "hbond_per_res",
    "disulfide_per_res",
    "sasa_per_res",
    "has_foldx_result",
    "foldx_source_layer",
]

for col in targets:
    if col not in df.columns:
        df[col] = np.nan

# Make this explicitly text/object so pandas accepts labels
df["foldx_source_layer"] = df["foldx_source_layer"].fillna("").astype(str)

# Force patch PDB rows using accession overlap
pdb_mask = df["structure_sources"].astype(str).str.strip().eq("pdb")
overlap_mask = pdb_mask & df["uniprot_accession"].isin(pdb_best.index)
patched = int(overlap_mask.sum())

for idx in df.index[overlap_mask]:
    acc = df.at[idx, "uniprot_accession"]
    r = pdb_best.loc[acc]

    df.at[idx, "foldx_structure_source"] = r.get("structure_source", "pdb")
    df.at[idx, "foldx_structure_id"] = r.get("structure_id", "")
    df.at[idx, "foldx_input_pdb"] = r.get("foldx_input_pdb", "")
    df.at[idx, "foldx_ready"] = r.get("foldx_ready", "")
    df.at[idx, "chain_length"] = r.get("effective_chain_length", np.nan)
    df.at[idx, "foldx_wt_total_energy"] = r.get("foldx_wt_total_energy", np.nan)
    df.at[idx, "foldx_energy_per_residue"] = r.get("foldx_energy_per_residue_calc", np.nan)
    df.at[idx, "hbond_proxy_count"] = r.get("hbond_count", np.nan)
    df.at[idx, "salt_bridge_count"] = r.get("salt_bridges", np.nan)
    df.at[idx, "disulfide_count"] = r.get("disulfides", np.nan)
    df.at[idx, "sasa_total"] = r.get("sasa", np.nan)
    df.at[idx, "hbond_per_res"] = r.get("hbond_per_res_calc", np.nan)
    df.at[idx, "disulfide_per_res"] = r.get("disulfide_per_res_calc", np.nan)
    df.at[idx, "sasa_per_res"] = r.get("sasa_per_res_calc", np.nan)
    df.at[idx, "foldx_source_layer"] = "pdb_foldx_with_structure"

# Recalculate FoldX flag
df["has_foldx_result"] = pd.to_numeric(df["foldx_wt_total_energy"], errors="coerce").notna()

# Restore modeller source labels
modeller_mask = df["structure_sources"].astype(str).str.strip().eq("modeller") & df["has_foldx_result"]
df.loc[modeller_mask, "foldx_source_layer"] = "modeller_foldx_normalized"

# Empty label for no-FoldX rows
df.loc[~df["has_foldx_result"], "foldx_source_layer"] = ""

df.to_csv(MASTER_OUT, index=False)

lines = []
lines.append("Final Xylanase Thesis Master Table v2 Audit")
lines.append("=" * 60)
lines.append(f"Input master: {MASTER_IN}")
lines.append(f"PDB FoldX source: {PDB_FOLDX}")
lines.append(f"Output: {MASTER_OUT}")
lines.append("")
lines.append(f"Rows: {len(df)}")
lines.append(f"Unique UniProt accessions: {df['uniprot_accession'].nunique()}")
lines.append(f"PDB rows in master: {int(pdb_mask.sum())}")
lines.append(f"PDB rows patched by accession overlap: {patched}")
lines.append("")
lines.append("Structure source summary:")
lines.append(str(df["structure_sources"].value_counts(dropna=False)))
lines.append("")
lines.append("FoldX result coverage:")
lines.append(str(df["has_foldx_result"].value_counts(dropna=False)))
lines.append("")
lines.append("FoldX source layer:")
lines.append(str(df["foldx_source_layer"].value_counts(dropna=False)))
lines.append("")
lines.append("FoldX coverage by structure source:")
lines.append(str(pd.crosstab(df["structure_sources"].fillna("none"), df["has_foldx_result"])))
lines.append("")
lines.append("Organism type x GH family:")
lines.append(str(pd.crosstab(df["organism_type"], df["gh_family"])))

AUDIT_OUT.write_text("\n".join(lines))

print(f"Saved: {MASTER_OUT}")
print(f"Saved: {AUDIT_OUT}")
print(f"Rows: {len(df)}")
print(f"Unique UniProt accessions: {df['uniprot_accession'].nunique()}")
print(f"PDB rows patched by accession overlap: {patched}")
print()
print(pd.crosstab(df["structure_sources"].fillna("none"), df["has_foldx_result"]))
print()
print(df["foldx_source_layer"].value_counts(dropna=False))
