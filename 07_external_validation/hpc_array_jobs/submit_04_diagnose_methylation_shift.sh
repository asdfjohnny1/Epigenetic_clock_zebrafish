#!/bin/bash
# SLURM submission script for 04_diagnose_methylation_shift.R - single job,
# not an array. Adjust #SBATCH resource directives for your own cluster.
# Requires hpc_array_jobs/submit_02_bismark_align_extract.sh to have
# completed already (reuses the same .cov.gz files as step 3).

#SBATCH --job-name=extval_diagnose_shift
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your-email@institution.ac.uk
#SBATCH -o log/extval_diagnose_shift-%j.out
#SBATCH -e log/extval_diagnose_shift-%j.err

module purge
module load gcc/11.1.0
module load R/4.5.1

srun Rscript "$EPICLOCK_ROOT/07_external_validation/04_diagnose_methylation_shift.R"
