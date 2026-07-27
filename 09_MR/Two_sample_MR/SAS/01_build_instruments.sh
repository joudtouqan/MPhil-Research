#!/bin/bash
#SBATCH --job-name=mr_instruments
#SBATCH --output=/home/jt962/rds/hpc-work/believe_adiposity/09_MR/logs/01_build_instruments_%j.out
#SBATCH --error=/home/jt962/rds/hpc-work/believe_adiposity/09_MR/logs/01_build_instruments_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH -p icelake
#SBATCH -A BUTTERWORTH-SL2-CPU

# ============================================================================
# Step 2 of the BELIEVE->G&H two-sample MR pipeline: build exposure instruments.
#
# Per trait (BMI, waist, hip, WHRadjBMI, fat mass, fat %):
#   1. Stream the 22 per-chrom SAIGE HDS sumstats, apply instrument filters:
#        - genome-wide significant   p.value < 5e-8
#        - common                    MAF (from AF_Allele2) >= 0.01
#        - well-imputed              imputationInfo >= 0.7
#        - exclude extended MHC      chr6:25,726,063-33,400,644 (GRCh38, xMHC)
#   2. LD-clump in PLINK2 against the in-sample BELIEVE reference 
#        --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.001 --clump-kb 10000
# ============================================================================

set -uo pipefail
module load rhel8/default-icl

PLINK2=/usr/local/Cluster-Apps/ceuadmin/plink/2.0_20240105/plink2
O=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/04_step2_saige/HDS/output
MR=/home/jt962/rds/hpc-work/believe_adiposity/09_MR
INST=$MR/instruments
BFILE=/home/jt962/rds/hpc-work/believe_adiposity/03_variant_QC/pgen/freeze_two_common_qc
MHC_LO=25726063
MHC_HI=33400644
mkdir -p "$INST"

# trait label -> SAIGE file stem
TRAITS="bmi:bmi_INT waist:waist_INT hip:hip_INT whradjbmi:whr_bmi_adj_INT fatmass:fatmass_INT fatperc:fatperc_INT"

SAIGE_HDR="CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\timputationInfo\tBETA\tSE\tTstat\tvar\tp.value\tN"

echo "=== Step 2: build instruments  $(date) ==="
SUMMARY=$INST/_instrument_counts.tsv
echo -e "trait\tn_gwsig_QC\tn_clumped_instruments" > "$SUMMARY"

for kv in $TRAITS; do
  name=${kv%%:*}; stem=${kv##*:}
  gw=$INST/${name}.gwsig.tsv
  ci=$INST/${name}.clumpin.tsv
  echo; echo "--- $name ($stem) ---"

  # headers (gwsig: SNP + the 14 SAIGE columns)
  printf "SNP\t%b\n" "$SAIGE_HDR" > "$gw"
  printf "ID\tP\n" > "$ci"

  # check all 22 chr files exist
  for c in $(seq 1 22); do
    [ -s "$O/${stem}_chr${c}.txt" ] || { echo "ERROR missing $O/${stem}_chr${c}.txt"; exit 1; }
  done

  # single awk pass over the 22 per-chrom files (FNR==1 skips each header)
  awk -F'\t' -v OFS='\t' -v lo=$MHC_LO -v hi=$MHC_HI -v gwf="$gw" -v cif="$ci" '
    FNR==1 { next }
    {
      af=$7+0; maf=(af<1-af)?af:1-af; info=$8+0; p=$13+0;
      if (p<5e-8 && maf>=0.01 && info>=0.7 && !($1==6 && $2>=lo && $2<=hi)) {
        id="chr"$3;
        print id, $0  >> gwf;
        print id, $13 >> cif;
      }
    }
  ' $O/${stem}_chr{1..22}.txt

  ngw=$(($(wc -l < "$gw") - 1))
  echo "  genome-wide-sig + QC variants (pre-clump): $ngw"

  if [ "$ngw" -eq 0 ]; then
    echo "  no variants pass filters; skipping clump"
    echo -e "${name}\t0\t0" >> "$SUMMARY"
    continue
  fi

  # LD clump against in-sample reference
  $PLINK2 --bfile "$BFILE" \
    --clump "$ci" \
    --clump-id-field ID --clump-p-field P \
    --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.001 --clump-kb 10000 \
    --out "$INST/${name}" \
    --threads 4 --memory 14000

  clumps=$INST/${name}.clumps
  if [ ! -s "$clumps" ]; then
    echo "  WARNING: no .clumps output for $name"
    echo -e "${name}\t${ngw}\tNA" >> "$SUMMARY"
    continue
  fi

  # index-variant IDs from .clumps (column named ID), then pull full rows from gwsig
  idcol=$(head -1 "$clumps" | tr '\t' '\n' | grep -nx "ID" | cut -d: -f1)
  awk -v c="$idcol" 'NR>1{print $c}' "$clumps" | sort -u > "$INST/${name}.index_ids.txt"
  nclump=$(wc -l < "$INST/${name}.index_ids.txt")

  inst=$INST/${name}.instruments.tsv
  head -1 "$gw" > "$inst"
  awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR==1{next} ($1 in keep)' \
      "$INST/${name}.index_ids.txt" "$gw" >> "$inst"
  ninst=$(($(wc -l < "$inst") - 1))

  echo "  clumped index variants (instruments): $ninst"
  echo -e "${name}\t${ngw}\t${ninst}" >> "$SUMMARY"
done

echo; echo "=== instrument counts ==="
cat "$SUMMARY"
echo; echo "=== done  $(date) ==="
