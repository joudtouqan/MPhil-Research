
set -euo pipefail

module load plink/2.00-alpha
module load singularity/current

# ── Paths ─────────────────────────────────────────────────────────────────────
LD=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/01_LD_pruning
WORK=/home/jt962/rds/hpc-work/believe_adiposity/05_GWAS/Saige/02_GRM
SIF=/home/jt962/rds/hpc-work/software/saige/saige_1.5.0.sif
mkdir -p $WORK/bfile $WORK/grm $WORK/logs

# ── Step A: convert LD-pruned PGEN → BED/BIM/FAM, autosomes only ─────────────
plink2 \
  --pfile $LD/ld_pruned_full_cohort \
  --autosome \
  --make-bed \
  --threads 20 \
  --memory 100000 \
  --out $WORK/bfile/ld_pruned_autosomes

# ── Step B: build sparse GRM via SAIGE container ─────────────────────────────
# Prevent the host ~/.Rprofile from leaking the cluster R libpath into the
# container, which would cause SAIGE to load the wrong (incompatible) Rcpp.
export APPTAINERENV_R_PROFILE_USER=/dev/null
export APPTAINERENV_R_ENVIRON_USER=/dev/null

singularity exec --cleanenv --bind /rds,/home $SIF \
  createSparseGRM.R \
  --plinkFile=$WORK/bfile/ld_pruned_autosomes \
  --nThreads=$SLURM_CPUS_PER_TASK \
  --outputPrefix=$WORK/grm/sparse_grm \
  --numRandomMarkerforSparseKin=2000 \
  --relatednessCutoff=0.05

echo "Sparse GRM build complete: $(date)"
ls -lh $WORK/grm/sparse_grm*

# ── Resource report ───────────────────────────────────────────────────────────
sleep 20
echo -e "\n####################################################" >&2
echo "RESOURCE USAGE REPORT for Job: $SLURM_JOB_ID" >&2
echo "Finished at: $(date)" >&2
echo "----------------------------------------------------" >&2
sacct -j $SLURM_JOB_ID --format=JobID,JobName%20,MaxRSS,TotalCPU,Elapsed,CPUTime,AllocCPUS,State >&2
echo "####################################################" >&2
