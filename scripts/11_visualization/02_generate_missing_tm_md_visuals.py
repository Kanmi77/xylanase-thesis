#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path.home() / "xylanase-thesis"

VIS_DIR = BASE / "results/visualization_thesis"
PRIORITY_DIR = VIS_DIR / "priority_13"
STRUCT_DIR = VIS_DIR / "structure_validation"
MD_DIR = VIS_DIR / "md"

for d in [VIS_DIR, PRIORITY_DIR, STRUCT_DIR, MD_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# 1. TM-align scatter: TM-score vs RMSD
# ---------------------------------------------------------------------
tm_file = BASE / "results/reports/tmalign_best_reference_quality_classes.csv"

if tm_file.exists():
    tm = pd.read_csv(tm_file)

    required = {"tm_score_best", "rmsd", "model_organism_type", "model_gh_family", "tmalign_quality_class"}
    missing = required - set(tm.columns)

    if missing:
        print(f"TM-align file found but missing columns: {missing}")
    else:
        tm["tm_score_best"] = pd.to_numeric(tm["tm_score_best"], errors="coerce")
        tm["rmsd"] = pd.to_numeric(tm["rmsd"], errors="coerce")
        tm = tm.dropna(subset=["tm_score_best", "rmsd"])

        plt.figure(figsize=(10, 7))

        for label, sub in tm.groupby("tmalign_quality_class"):
            plt.scatter(
                sub["rmsd"],
                sub["tm_score_best"],
                s=28,
                alpha=0.70,
                label=f"{label} (n={len(sub)})"
            )

        plt.axhline(0.95, linestyle="--", linewidth=1)
        plt.axvline(1.0, linestyle="--", linewidth=1)

        plt.xlabel("RMSD to best structural reference")
        plt.ylabel("Best TM-score")
        plt.title("TM-align structural validation of modelled and experimental xylanase structures")
        plt.legend(frameon=False, fontsize=9)
        plt.tight_layout()

        out1 = STRUCT_DIR / "tm_score_vs_rmsd.png"
        out2 = PRIORITY_DIR / "03_tm_score_vs_rmsd.png"

        plt.savefig(out1, dpi=300)
        plt.savefig(out2, dpi=300)
        plt.close()

        print(f"Saved: {out1}")
        print(f"Saved: {out2}")

        # Extra grouped bar: quality class distribution
        quality_counts = (
            tm.groupby(["model_organism_type", "model_gh_family", "tmalign_quality_class"])
            .size()
            .reset_index(name="count")
        )
        quality_counts.to_csv(STRUCT_DIR / "tmalign_quality_distribution_for_plot.csv", index=False)

else:
    print(f"TM-align source file missing: {tm_file}")


# ---------------------------------------------------------------------
# 2. MD status visual from md_tier2_wt_mutant_compact
# ---------------------------------------------------------------------
md_base = BASE / "md_tier2_wt_mutant_compact/systems"

rows = []

if md_base.exists():
    for system_dir in sorted(md_base.iterdir()):
        if not system_dir.is_dir():
            continue

        for temp_dir in sorted(system_dir.iterdir()):
            if not temp_dir.is_dir():
                continue

            temp = temp_dir.name
            status_file = temp_dir / "status.txt"
            md_log = temp_dir / "md.log"
            md_xtc = temp_dir / "md.xtc"
            md_gro = temp_dir / "md.gro"
            md_cpt = temp_dir / "md.cpt"

            status_text = ""
            if status_file.exists():
                try:
                    status_text = status_file.read_text(errors="ignore").strip()
                except Exception:
                    status_text = ""

            has_log = md_log.exists()
            has_xtc = md_xtc.exists()
            has_gro = md_gro.exists()
            has_cpt = md_cpt.exists()

            if "DONE" in status_text.upper() or "COMPLETE" in status_text.upper():
                status_numeric = 3
                status_label = "completed"
            elif has_log and has_xtc and has_cpt:
                status_numeric = 2
                status_label = "running_or_interrupted_with_outputs"
            elif has_log or has_xtc or has_cpt:
                status_numeric = 1
                status_label = "started_partial_outputs"
            else:
                status_numeric = 0
                status_label = "not_started"

            rows.append({
                "system": system_dir.name,
                "temperature": temp,
                "status": status_label,
                "status_numeric": status_numeric,
                "has_md_log": has_log,
                "has_md_xtc": has_xtc,
                "has_md_gro": has_gro,
                "has_md_cpt": has_cpt,
                "status_text": status_text,
            })

    md_status = pd.DataFrame(rows)
    md_status_file = MD_DIR / "md_tier2_status_table.csv"
    md_status.to_csv(md_status_file, index=False)
    print(f"Saved: {md_status_file}")

    if not md_status.empty:
        pivot = md_status.pivot_table(
            index="system",
            columns="temperature",
            values="status_numeric",
            aggfunc="max",
            fill_value=0
        )

        plt.figure(figsize=(11, max(5, 0.42 * len(pivot))))
        plt.imshow(pivot.values, aspect="auto")

        plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=0)
        plt.yticks(range(len(pivot.index)), pivot.index, fontsize=8)

        cbar = plt.colorbar()
        cbar.set_ticks([0, 1, 2, 3])
        cbar.set_ticklabels([
            "not started",
            "partial",
            "outputs present",
            "completed"
        ])

        plt.xlabel("Simulation temperature")
        plt.ylabel("WT/mutant MD system")
        plt.title("Current MD simulation status for WT–mutant validation systems")
        plt.tight_layout()

        out3 = MD_DIR / "md_tier2_simulation_status_heatmap.png"
        out4 = PRIORITY_DIR / "12_md_simulation_status_heatmap.png"

        plt.savefig(out3, dpi=300)
        plt.savefig(out4, dpi=300)
        plt.close()

        print(f"Saved: {out3}")
        print(f"Saved: {out4}")

else:
    print(f"MD system directory missing: {md_base}")


# ---------------------------------------------------------------------
# 3. Final priority directory audit
# ---------------------------------------------------------------------
priority_files = sorted(PRIORITY_DIR.glob("*.png"))

audit = PRIORITY_DIR / "priority_13_visual_audit.txt"
with audit.open("w") as f:
    f.write("Priority thesis visual audit\n")
    f.write("=" * 60 + "\n")
    f.write(f"Directory: {PRIORITY_DIR}\n")
    f.write(f"PNG files found: {len(priority_files)}\n\n")
    for p in priority_files:
        f.write(f"- {p.name}\n")

print(f"Saved: {audit}")
print("")
print("Priority files now present:")
for p in priority_files:
    print(p.name)
