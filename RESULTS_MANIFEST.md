# Results Manifest for Thesis Analysis

This file lists the result files selected for GitHub tracking and thesis analysis.

## 1. Data Acquisition and Curation

| File Path | Description |
|---|---|
| `data/curated/xylanase_master_all_curated_with_brenda.csv` | Final master table with GH families, organism types, and BRENDA temperature/pH annotations. |
| `data/curated/xylanase_master_deduplicated.csv` | Deduplicated master table, one entry per UniProt accession. |
| `data/curated/refseq_inventory.csv` | Mapping of RefSeq accessions to UniProt, organism, and sequence. |
| `data/curated/foldx_combined_inputs.csv` | Combined list of all structures with paths to FoldX-ready PDBs. |

## 2. Sequence and Phylogenetic Analysis

| File Path | Description |
|---|---|
| `results/alignments/all_GH10_GH11.mafft.fasta` | MAFFT alignment of all GH10/GH11 sequences. |
| `results/alignments/bacterial_GH10.mafft.fasta` | MAFFT alignment for bacterial GH10. |
| `results/alignments/bacterial_GH11.mafft.fasta` | MAFFT alignment for bacterial GH11. |
| `results/alignments/fungal_GH10.mafft.fasta` | MAFFT alignment for fungal GH10. |
| `results/alignments/fungal_GH11.mafft.fasta` | MAFFT alignment for fungal GH11. |
| `results/alignments/archaea_GH10.mafft.fasta` | MAFFT alignment for archaea GH10. |
| `results/alignments/archaea_GH11.mafft.fasta` | MAFFT alignment for archaea GH11. |
| `results/trees/all_GH10_GH11.fasttree.nwk` | Phylogenetic tree for combined GH10/GH11 dataset. |
| `results/trees/bacterial_GH10.fasttree.nwk` | Phylogenetic tree for bacterial GH10. |
| `results/trees/bacterial_GH11.fasttree.nwk` | Phylogenetic tree for bacterial GH11. |
| `results/trees/fungal_GH10.fasttree.nwk` | Phylogenetic tree for fungal GH10. |
| `results/trees/fungal_GH11.fasttree.nwk` | Phylogenetic tree for fungal GH11. |
| `results/sequence_features/protparam_features.csv` | Physicochemical features for curated proteins. |
| `results/sequence_features/signal_peptide_proxy.csv` | Secretion signal proxy for curated proteins. |
| `results/sequence_features/consensus/bacterial_GH10_consensus.csv` | Conservation scores for bacterial GH10 alignment. |
| `results/sequence_features/consensus/bacterial_GH11_consensus.csv` | Conservation scores for bacterial GH11 alignment. |
| `results/sequence_features/consensus/fungal_GH10_consensus.csv` | Conservation scores for fungal GH10 alignment. |
| `results/sequence_features/consensus/fungal_GH11_consensus.csv` | Conservation scores for fungal GH11 alignment. |
| `results/sequence_features/consensus/consensus_summary.csv` | Summary of conserved positions per group. |
| `results/sequence_features/consensus_stabilizing_candidates.csv` | Conserved residues merged across groups. |
| `results/sequence_features/GH10_vs_GH11_signature_positions.csv` | Alignment positions where GH10 and GH11 differ in dominant residue. |
| `results/sequence_features/candidate_motif_positions.csv` | Putative motif positions mapped to alignments. |

## 3. Structural Characterization and Homology Modeling

| File Path | Description |
|---|---|
| `results/structures/combined_structural_features.csv` | Structural features for experimental and modelled structures. |
| `results/structures/tmalign_results_GH10_GH11_only.csv` | TM-align summary scores for model-reference pairs. |
| `results/structures/tmalign_best_reference_per_model.csv` | Best reference structure for each model. |
| `results/structures/tmalign_results_filtered.csv` | Filtered TM-align pairs used for residue mapping. |
| `results/foldx/foldx_normalized.csv` | Normalized FoldX energy per residue with metadata. |
| `data/curated/combined_structure_manifest.csv` | Master list of all structures with paths and sources. |
| `data/curated/modeller_model_manifest.csv` | MODELLER model manifest with template and quality information. |

## 4. Thermostability and Mutation Prediction

| File Path | Description |
|---|---|
| `results/foldx_clean/tier2_ddg_ranked.csv` | Raw ΔΔG values for Tier-2 mutations from FoldX BuildModel. |
| `results/foldx_clean/tier2_ddg_ranked_annotated.csv` | ΔΔG table annotated with mutation string, metadata, and industrial properties. |
| `results/foldx_clean/tier2_ddg_summary_by_protein.csv` | Per-protein mutation summary. |
| `results/foldx_clean/tier2_ddg_summary_by_group.csv` | Mutation summary by GH family and organism type. |
| `results/foldx_clean/tier2_top_stabilizing_mutations.csv` | Top stabilising mutations sorted by ΔΔG. |
| `data/curated/foldx_mutation_lists_tier2_fixed/` | FoldX-valid Tier-2 mutation files. |
| `results/ml/structural_stability_ml_predictions.csv` | Random forest cross-validation predictions. |
| `results/ml/structural_stability_ml_summary.csv` | ML model performance summary. |

## 5. Molecular Dynamics Validation

| File Path | Description |
|---|---|
| `md_tier2_wt_mutant/inputs/` | Wild-type and mutant PDB files staged for MD simulations. |

Large MD trajectory/output files are intentionally excluded from this manifest unless specifically selected later.

## 6. Workflow Output and Integration

| File Path | Description |
|---|---|
| `results/foldx/foldx_structural_ranked.csv` | Structures ranked by composite structural stability score. |
| `results/foldx_clean/tier2_ddg_ranked_industrial_annotated.csv` | Final integrated Tier-2 mutation dataset with BRENDA and secretion proxy annotations. |
