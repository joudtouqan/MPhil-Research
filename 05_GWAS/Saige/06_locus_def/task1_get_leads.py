#!/usr/bin/env python
"""TASK 1 - Distance-based locus definition (lead SNPs) for the 6 BELIEVE HDS
adiposity traits, on the FULL UNFILTERED HDS sumstats.

get_lead(sig_level=5e-8, windowsizekb=500): a new locus starts when the gap to
the previous significant variant exceeds 500 kb; the lead is the lowest-P SNP in
each locus. This defines locus boundaries/leads only -- NOT independence
(handled later by conditional analysis).

EA = Allele2 (SAIGE reports BETA & AF w.r.t. Allele2 regardless of --AlleleOrder).
"""
import os, sys
import pandas as pd
import gwaslab as gl

DEST = "/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/05_post_gwas/HDS"
SUM  = os.path.join(DEST, "sumstats")
TRAITS = ["bmi_INT","fatmass_INT","fatperc_INT","hip_INT",
          "waist_INT","whr_bmi_adj_INT"]
WIN_KB = 500
HALF   = WIN_KB * 1000 // 2     # nominal +/- flank for window_start/end

print("gwaslab", gl.__version__, flush=True)
summary = []
all_leads = []
for t in TRAITS:
    print(f"\n===== {t} =====", flush=True)
    ss = gl.Sumstats(
        os.path.join(SUM, f"{t}.sumstats.tsv"),
        sep="\t",
        snpid="SNP", chrom="CHR", pos="POS",
        ea="Allele2", nea="Allele1",          # EA = Allele2
        eaf="AF", beta="BETA", se="SE", p="P", n="N",
        other=["INFO","MAF"], build="38", verbose=True)

    leads = ss.get_lead(sig_level=5e-8, windowsizekb=WIN_KB,
                        verbose=True).copy()
    leads["TRAIT"] = t
    leads["WINDOW_START"] = (leads["POS"] - HALF).clip(lower=1)
    leads["WINDOW_END"]   = leads["POS"] + HALF
    out = os.path.join(DEST, f"{t}.leads.tsv")
    leads.to_csv(out, sep="\t", index=False)
    n = len(leads)
    print(f"{t}: {n} leads -> {out}", flush=True)
    summary.append((t, n))
    all_leads.append(leads)

# combined leads table
allleads = pd.concat(all_leads, ignore_index=True)
allleads.to_csv(os.path.join(DEST, "all_traits.leads.tsv"), sep="\t", index=False)

print("\n===== LEAD COUNT PER TRAIT =====", flush=True)
for t, n in summary:
    print(f"  {t:18s} {n}", flush=True)

# FTO sanity (chr16 ~53.7-54.2 Mb)
print("\n===== FTO sanity (chr16 53.7-54.2 Mb leads) =====", flush=True)
fto = allleads[(allleads.CHR==16) & (allleads.POS>=53_700_000)
               & (allleads.POS<=54_200_000)]
cols = [c for c in ["TRAIT","CHR","POS","SNPID","EA","NEA","EAF","P"] if c in fto.columns]
print(fto[cols].to_string(index=False) if len(fto) else "  (no FTO-window leads)",
      flush=True)
