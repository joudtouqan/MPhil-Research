#!/usr/bin/env python3
"""Summarise the HDS FINEMAP: iterate every (trait, lead) window in manifest.tsv,
reusing the original parse_finemap.py credible-set / PIP / log10-BF parsing so numbers stay
comparable to the six-locus run. Writes results/SUMMARY_sweep_hds.tsv.

Lead PIP is matched on the exact variant id (z-file 'rsid' == chr-prefixed SNPID) rather than on
chromosome:position. 
"""
import os, csv

D = "/home/jt962/rds/hpc-work/believe_adiposity/11_finemap/sweep_hds"
MANIFEST = f"{D}/manifest.tsv"

# ---- parsing helpers: identical logic to parse_finemap.py (lead matched by id) ----
def lead_pip(snp_file, lead_id):
    if not os.path.exists(snp_file): return (None, None, None)
    with open(snp_file) as f:
        r = csv.DictReader(f, delimiter=" ")
        rows = [x for x in r if x.get("rsid")]
    rows = [x for x in rows if x["rsid"]]
    lp = None
    for x in rows:
        if x["rsid"] == lead_id:
            lp = float(x["prob"]); break
    top = sorted(rows, key=lambda x: -float(x["prob"]))[:3]
    return lp, [(x["rsid"], round(float(x["prob"]),3), round(float(x["z"]),1)) for x in top], len(rows)

def cred_sizes(cred_file):
    if not os.path.exists(cred_file): return []
    with open(cred_file) as f:
        hdr = f.readline().split()
        sets = [h for h in hdr if h.startswith("cred_set")]
        n = len(sets)
        rows = [ln.split() for ln in f if ln.strip()]
    sizes = []
    for si in range(n):
        col = 1 + si*2
        cnt = sum(1 for rr in rows if len(rr)>col and rr[col] not in ("NA",""))
        sizes.append(cnt)
    return sizes

def logbf(log_file):
    if not os.path.exists(log_file): return None
    for ln in open(log_file):
        if "Log10-BF of >= one causal" in ln:
            return ln.split(":")[-1].strip()
    return None

with open(MANIFEST) as f:
    man = list(csv.DictReader(f, delimiter="\t"))

rows_out = []
low_info = []       # leads with INFO < 0.6  (flag only)
lead_absent = []    # leads not present in the LD reference (indels)
n_done = 0; n_missing = 0
print(f'{"tag":<22}{"trait":<16}{"k":>2}{"ncred":>6}{"cred_sizes":>14}{"leadPIP":>9}{"INFO":>7}  log10BF')
print("-"*112)
for m in man:
    tag, trait, k = m["tag"], m["trait"], m["ncausal"]
    lead_id = m["lead_snp"]
    try: info = float(m["lead_info"])
    except (TypeError, ValueError): info = None
    info_str = f"{info:.3f}" if info is not None else "NA"
    if info is not None and info < 0.6:
        low_info.append((tag, trait, lead_id, info))

    snp = f"{D}/results/{tag}.snp"; cred = f"{D}/results/{tag}.cred"; log = f"{D}/results/{tag}.log_sss"
    if not os.path.exists(log):
        n_missing += 1
        rows_out.append([tag, trait, lead_id, k, "SKIPPED/NO_OUTPUT", "", "", "", "", info_str, ""])
        continue
    n_done += 1
    lp, top, nv = lead_pip(snp, lead_id)
    sizes = cred_sizes(cred); bf = logbf(log)
    if lp is None:
        lead_absent.append((tag, trait, lead_id))
    print(f'{tag:<22}{trait:<16}{k:>2}{len(sizes):>6}{str(sizes):>14}'
          f'{(round(lp,3) if lp is not None else "NA"):>9}{info_str:>7}  {bf}')
    rows_out.append([tag, trait, lead_id, k, len(sizes), sizes,
                     (round(lp,4) if lp is not None else "NA"), bf, nv, info_str, top])

out_tsv = f"{D}/results/SUMMARY_sweep_hds.tsv"
with open(out_tsv, "w") as o:
    o.write("tag\ttrait\tlead_variant\tk_ncausal\tn_credible_sets\tcred_set_sizes\t"
            "lead_PIP\tlog10BF_ge1causal\tn_variants_window\tlead_INFO\ttop3_variants(id,PIP,z)\n")
    for r in rows_out:
        o.write("\t".join(str(x) for x in r)+"\n")

print("\n" + "="*70)
print(f"windows in manifest : {len(man)}")
print(f"windows with output : {n_done}")
print(f"windows missing out : {n_missing}")

print(f"\nLeads with imputation INFO < 0.6 (FLAGGED, not dropped): {len(low_info)}")
for tag, trait, snp_, info in sorted(low_info, key=lambda x: x[3]):
    print(f"  {tag:<22} {trait:<16} {snp_:<26} INFO={info:.3f}")

print(f"\nLeads absent from the SNP-only LD reference (lead PIP = NA, window still fine-mapped): {len(lead_absent)}")
for tag, trait, snp_ in lead_absent:
    print(f"  {tag:<22} {trait:<16} {snp_}")

print(f"\nWrote {out_tsv}")
