#!/bin/bash
# SLURM submission script for 03_apply_best_clocks.R - single job, not an
# array (reads all external samples' coverage files together). Adjust
# #SBATCH resource directives for your own cluster. Requires
# hpc_array_jobs/submit_02_bismark_align_extract.sh to have completed for
# every sample listed in external_metadata_template.csv first.

#SBATCH --job-name=extval_apply_clocks
#SBATCH --partition=compute
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your-email@institution.ac.uk
#SBATCH -o log/extval_apply_clocks-%j.out
#SBATCH -e log/extval_apply_clocks-%j.err

module purge
module load gcc/11.1.0
module load R/4.5.1

srun Rscript "$EPICLOCK_ROOT/07_external_validation/03_apply_best_clocks.R"
