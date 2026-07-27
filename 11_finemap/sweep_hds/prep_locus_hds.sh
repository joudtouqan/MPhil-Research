#!/bin/bash
# Prepare FINEMAP inputs for one locus with in-sample BELIEVE LD -- HDS sumstats variant.
# Usage: prep_locus_hds.sh <trait> <chr> <index_pos> <tag>
# Produces z/<tag>.z, ld/<tag>.ld (same variant order, effect allele = ALT = SAIGE Allele2).

set -euo pipefail
trait=$1; chr=$2; pos=$3; tag=$4
WIN=250000
ROOT=/home/jt962/rds/hpc-work/believe_adiposity
PV=$ROOT/03_variant_QC/pgen/freeze_two_common_qc           # bed/bim/fam
SS=$ROOT/05_GWAS/Saige/05_post_gwas/HDS/sumstats/${trait}.sumstats.tsv
PLINK=/usr/local/Cluster-Apps/plink/1.9/plink
OUT=$ROOT/11_finemap/sweep_hds
lo=$((pos-WIN)); hi=$((pos+WIN)); [ $lo -lt 1 ] && lo=1

# 1. candidate variants in window from sumstats with MAF>=0.01; prefix id to match pvar/bim
awk -F'\t' -v c=$chr -v lo=$lo -v hi=$hi \
  'NR>1 && $1==c && $2>=lo && $2<=hi && $14>=0.01 {print "chr"$3}' "$SS" \
  | sort -u > $OUT/z/${tag}.cand.txt
nc=$(wc -l < $OUT/z/${tag}.cand.txt)
echo "[$tag] candidate variants in +-${WIN}bp window: $nc"

# 2. retained variant order from plink (keep-allele-order: A1 stays = ALT = effect allele)
$PLINK --bfile $PV --extract $OUT/z/${tag}.cand.txt --keep-allele-order \
       --make-just-bim --out $OUT/z/${tag} >/dev/null 2>&1
nr=$(wc -l < $OUT/z/${tag}.bim)
echo "[$tag] variants present in reference: $nr"

# 3. signed correlation square matrix (space-delimited) in .bim order
$PLINK --bfile $PV --extract $OUT/z/${tag}.cand.txt --keep-allele-order \
       --r square spaces --out $OUT/ld/${tag} >/dev/null 2>&1

# 4. build z-file in EXACT .bim order; effect allele = bim A1 (col5) = ALT
#    join sumstat beta/se/maf by chr-prefixed variant id
awk -F'\t' -v ss="$SS" '
  BEGIN{ OFS=" ";
        while((getline line < ss)>0){ n=split(line,a,"\t");
          if(a[1]=="CHR") continue; id="chr"a[3];
          b[id]=a[9]; se[id]=a[10]; maf[id]=a[14]; }
        print "rsid","chromosome","position","allele1","allele2","maf","beta","se"; }
  { id=$2; chr=$1; bp=$4; a1=$5; a2=$6;   # bim: a1=ALT(effect), a2=REF
    if(id in b){ print id, chr, bp, a1, a2, maf[id], b[id], se[id] } }
' $OUT/z/${tag}.bim > $OUT/z/${tag}.z
nz=$(($(wc -l < $OUT/z/${tag}.z)-1))
echo "[$tag] z-file variants: $nz   (z rows must equal LD dim)"
echo "[$tag] LD matrix dim: $(head -1 $OUT/ld/${tag}.ld | awk '{print NF}') x $(wc -l < $OUT/ld/${tag}.ld)"
