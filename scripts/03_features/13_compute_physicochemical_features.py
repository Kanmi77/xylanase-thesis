#!/usr/bin/env python3
"""
Compute ProtParam-like physicochemical features from sequences.

Input:
  data/curated/xylanase_master_all_curated.csv

Output:
  results/features/physicochemical_features.csv

Features:
- length
- aa composition (fractions)
- molecular weight (Da)
- theoretical pI (approx; Bjellqvist pK set)
- GRAVY (Kyte-Doolittle)
- aromaticity (F+W+Y)/length
- instability index (Guruprasad et al. dipeptide weights)
- aliphatic index (Ikai)
- counts: acidic (D+E), basic (K+R+H)
"""

import os, math
import pandas as pd

BASE = os.path.expanduser("~/xylanase-thesis")
IN_CSV = os.path.join(BASE, "data/curated/xylanase_master_all_curated.csv")
OUT_CSV = os.path.join(BASE, "results/features/physicochemical_features.csv")

AA = list("ACDEFGHIKLMNPQRSTVWY")

# Average residue masses (monoisotopic is overkill; use average residue weights)
RES_MASS = {
    "A": 71.0788, "C": 103.1388, "D": 115.0886, "E": 129.1155, "F": 147.1766,
    "G": 57.0519, "H": 137.1411, "I": 113.1594, "K": 128.1741, "L": 113.1594,
    "M": 131.1926, "N": 114.1038, "P": 97.1167, "Q": 128.1307, "R": 156.1875,
    "S": 87.0782, "T": 101.1051, "V": 99.1326, "W": 186.2132, "Y": 163.1760
}
WATER = 18.01528

# Kyte-Doolittle hydropathy
KD = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3
}

# pK values (Bjellqvist-like)
PK_NTERM = 9.69
PK_CTERM = 2.34
PK_SIDE = {"C": 8.33, "D": 3.86, "E": 4.25, "H": 6.00, "K": 10.53, "R": 12.48, "Y": 10.07}

# Instability index dipeptide weights (Guruprasad et al.)
# Full 400-table is large; we store it externally? No — we include as a compact dict-of-dicts.
# Source values are standard ProtParam table. This is the minimal reproducible way.
# NOTE: This table is long but deterministic; keep it as provided here.
II = {
"AA":1.0,"AC":1.0,"AD":-7.0,"AE":1.0,"AF":1.0,"AG":-7.0,"AH":1.0,"AI":1.0,"AK":1.0,"AL":1.0,"AM":1.0,"AN":1.0,"AP":1.0,"AQ":1.0,"AR":1.0,"AS":1.0,"AT":1.0,"AV":1.0,"AW":1.0,"AY":1.0,
"CA":1.0,"CC":1.0,"CD":20.0,"CE":1.0,"CF":1.0,"CG":1.0,"CH":33.0,"CI":1.0,"CK":1.0,"CL":1.0,"CM":1.0,"CN":1.0,"CP":20.0,"CQ":1.0,"CR":1.0,"CS":1.0,"CT":1.0,"CV":1.0,"CW":1.0,"CY":1.0,
"DA":1.0,"DC":1.0,"DD":1.0,"DE":1.0,"DF":1.0,"DG":1.0,"DH":1.0,"DI":1.0,"DK":1.0,"DL":1.0,"DM":1.0,"DN":1.0,"DP":1.0,"DQ":1.0,"DR":1.0,"DS":1.0,"DT":1.0,"DV":1.0,"DW":1.0,"DY":1.0,
"EA":1.0,"EC":1.0,"ED":1.0,"EE":1.0,"EF":1.0,"EG":1.0,"EH":1.0,"EI":1.0,"EK":1.0,"EL":1.0,"EM":1.0,"EN":1.0,"EP":1.0,"EQ":1.0,"ER":1.0,"ES":1.0,"ET":1.0,"EV":1.0,"EW":1.0,"EY":1.0,
"FA":1.0,"FC":1.0,"FD":13.0,"FE":1.0,"FF":1.0,"FG":1.0,"FH":1.0,"FI":1.0,"FK":1.0,"FL":1.0,"FM":1.0,"FN":1.0,"FP":20.0,"FQ":1.0,"FR":1.0,"FS":1.0,"FT":1.0,"FV":1.0,"FW":1.0,"FY":1.0,
"GA":1.0,"GC":1.0,"GD":1.0,"GE":1.0,"GF":1.0,"GG":1.0,"GH":1.0,"GI":1.0,"GK":1.0,"GL":1.0,"GM":1.0,"GN":1.0,"GP":1.0,"GQ":1.0,"GR":1.0,"GS":1.0,"GT":1.0,"GV":1.0,"GW":1.0,"GY":1.0,
"HA":1.0,"HC":1.0,"HD":1.0,"HE":1.0,"HF":1.0,"HG":1.0,"HH":1.0,"HI":1.0,"HK":1.0,"HL":1.0,"HM":1.0,"HN":1.0,"HP":1.0,"HQ":1.0,"HR":1.0,"HS":1.0,"HT":1.0,"HV":1.0,"HW":1.0,"HY":1.0,
"IA":1.0,"IC":1.0,"ID":1.0,"IE":1.0,"IF":1.0,"IG":1.0,"IH":1.0,"II":1.0,"IK":1.0,"IL":1.0,"IM":1.0,"IN":1.0,"IP":1.0,"IQ":1.0,"IR":1.0,"IS":1.0,"IT":1.0,"IV":1.0,"IW":1.0,"IY":1.0,
"KA":1.0,"KC":1.0,"KD":1.0,"KE":1.0,"KF":1.0,"KG":1.0,"KH":1.0,"KI":1.0,"KK":1.0,"KL":1.0,"KM":1.0,"KN":1.0,"KP":1.0,"KQ":1.0,"KR":1.0,"KS":1.0,"KT":1.0,"KV":1.0,"KW":1.0,"KY":1.0,
"LA":1.0,"LC":1.0,"LD":1.0,"LE":1.0,"LF":1.0,"LG":1.0,"LH":1.0,"LI":1.0,"LK":1.0,"LL":1.0,"LM":1.0,"LN":1.0,"LP":1.0,"LQ":1.0,"LR":1.0,"LS":1.0,"LT":1.0,"LV":1.0,"LW":1.0,"LY":1.0,
"MA":1.0,"MC":1.0,"MD":1.0,"ME":1.0,"MF":1.0,"MG":1.0,"MH":1.0,"MI":1.0,"MK":1.0,"ML":1.0,"MM":1.0,"MN":1.0,"MP":1.0,"MQ":1.0,"MR":1.0,"MS":1.0,"MT":1.0,"MV":1.0,"MW":1.0,"MY":1.0,
"NA":1.0,"NC":1.0,"ND":1.0,"NE":1.0,"NF":1.0,"NG":1.0,"NH":1.0,"NI":1.0,"NK":1.0,"NL":1.0,"NM":1.0,"NN":1.0,"NP":1.0,"NQ":1.0,"NR":1.0,"NS":1.0,"NT":1.0,"NV":1.0,"NW":1.0,"NY":1.0,
"PA":1.0,"PC":1.0,"PD":1.0,"PE":1.0,"PF":1.0,"PG":1.0,"PH":1.0,"PI":1.0,"PK":1.0,"PL":1.0,"PM":1.0,"PN":1.0,"PP":1.0,"PQ":1.0,"PR":1.0,"PS":1.0,"PT":1.0,"PV":1.0,"PW":1.0,"PY":1.0,
"QA":1.0,"QC":1.0,"QD":1.0,"QE":1.0,"QF":1.0,"QG":1.0,"QH":1.0,"QI":1.0,"QK":1.0,"QL":1.0,"QM":1.0,"QN":1.0,"QP":1.0,"QQ":1.0,"QR":1.0,"QS":1.0,"QT":1.0,"QV":1.0,"QW":1.0,"QY":1.0,
"RA":1.0,"RC":1.0,"RD":1.0,"RE":1.0,"RF":1.0,"RG":1.0,"RH":1.0,"RI":1.0,"RK":1.0,"RL":1.0,"RM":1.0,"RN":1.0,"RP":1.0,"RQ":1.0,"RR":1.0,"RS":1.0,"RT":1.0,"RV":1.0,"RW":1.0,"RY":1.0,
"SA":1.0,"SC":1.0,"SD":1.0,"SE":1.0,"SF":1.0,"SG":1.0,"SH":1.0,"SI":1.0,"SK":1.0,"SL":1.0,"SM":1.0,"SN":1.0,"SP":1.0,"SQ":1.0,"SR":1.0,"SS":1.0,"ST":1.0,"SV":1.0,"SW":1.0,"SY":1.0,
"TA":1.0,"TC":1.0,"TD":1.0,"TE":1.0,"TF":1.0,"TG":1.0,"TH":1.0,"TI":1.0,"TK":1.0,"TL":1.0,"TM":1.0,"TN":1.0,"TP":1.0,"TQ":1.0,"TR":1.0,"TS":1.0,"TT":1.0,"TV":1.0,"TW":1.0,"TY":1.0,
"VA":1.0,"VC":1.0,"VD":1.0,"VE":1.0,"VF":1.0,"VG":1.0,"VH":1.0,"VI":1.0,"VK":1.0,"VL":1.0,"VM":1.0,"VN":1.0,"VP":1.0,"VQ":1.0,"VR":1.0,"VS":1.0,"VT":1.0,"VV":1.0,"VW":1.0,"VY":1.0,
"WA":1.0,"WC":1.0,"WD":1.0,"WE":1.0,"WF":1.0,"WG":1.0,"WH":1.0,"WI":1.0,"WK":1.0,"WL":1.0,"WM":1.0,"WN":1.0,"WP":1.0,"WQ":1.0,"WR":1.0,"WS":1.0,"WT":1.0,"WV":1.0,"WW":1.0,"WY":1.0,
"YA":1.0,"YC":1.0,"YD":1.0,"YE":1.0,"YF":1.0,"YG":1.0,"YH":1.0,"YI":1.0,"YK":1.0,"YL":1.0,"YM":1.0,"YN":1.0,"YP":1.0,"YQ":1.0,"YR":1.0,"YS":1.0,"YT":1.0,"YV":1.0,"YW":1.0,"YY":1.0
}

# Note: The real ProtParam instability index uses a specific 400-value table.
# Here we've set non-provided pairs to 1.0 by default, which is NOT identical.
# To stay faithful, we implement the standard table via BioPython below instead.

from Bio.SeqUtils.ProtParam import ProteinAnalysis

def compute(seq: str) -> dict:
    seq = seq.strip().upper()
    seq = "".join([aa for aa in seq if aa in AA])
    if not seq:
        return {}

    pa = ProteinAnalysis(seq)

    length = len(seq)
    comp = pa.count_amino_acids()
    comp_frac = {f"aa_frac_{a}": comp.get(a,0)/length for a in AA}

    # ProtParam-like metrics from BioPython (closest reproducible analog)
    mw = pa.molecular_weight()
    gravy = pa.gravy()
    arom = pa.aromaticity()
    instab = pa.instability_index()
    try:
        pi = pa.isoelectric_point()
    except Exception:
        pi = float("nan")

    # Aliphatic index (Ikai)
    # AI = X(Ala)*100 + a*X(Val)*100 + b*(X(Ile)+X(Leu))*100
    # a=2.9, b=3.9
    xA = comp.get("A",0)/length
    xV = comp.get("V",0)/length
    xI = comp.get("I",0)/length
    xL = comp.get("L",0)/length
    aliphatic_index = (xA + 2.9*xV + 3.9*(xI + xL)) * 100.0

    acidic = comp.get("D",0) + comp.get("E",0)
    basic = comp.get("K",0) + comp.get("R",0) + comp.get("H",0)

    out = {
        "length": length,
        "molecular_weight": mw,
        "theoretical_pI": pi,
        "gravy": gravy,
        "aromaticity": arom,
        "instability_index": instab,
        "aliphatic_index": aliphatic_index,
        "acidic_count_DE": acidic,
        "basic_count_KRH": basic,
    }
    out.update(comp_frac)
    return out

def main() -> int:
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df = pd.read_csv(IN_CSV, dtype=str).fillna("")

    feats = []
    bad = 0

    for _, r in df.iterrows():
        seq = r.get("sequence","")
        f = compute(seq)
        if not f:
            bad += 1
            continue
        f.update({
            "primary_id": r.get("primary_id",""),
            "source": r.get("source",""),
            "uniprot_accession": r.get("uniprot_accession",""),
            "refseq_acc": r.get("refseq_acc",""),
            "organism_type": r.get("organism_type",""),
            "organism": r.get("organism",""),
            "gh_family": r.get("gh_family",""),
        })
        feats.append(f)

    out = pd.DataFrame(feats)
    out.to_csv(OUT_CSV, index=False)

    print(f"Wrote: {OUT_CSV} rows={len(out)} bad_sequences_skipped={bad}")
    print("\nFeature columns:", len(out.columns))
    print("\nCounts by organism_type:")
    print(out["organism_type"].value_counts().to_string())
    print("\nCounts by gh_family:")
    print(out["gh_family"].value_counts().to_string())

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
