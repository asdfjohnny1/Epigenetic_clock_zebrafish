#!/bin/bash
# SLURM array submission script - launches one task per sample listed in
# 07_external_validation/sample_list.txt. Adjust #SBATCH resource
# directives for your own cluster.

#SBATCH --job-name=extval_trim
#SBATCH --partition=compute
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your-email@institution.ac.uk
#SBATCH -o log/extval_trim-%A-%a.out
#SBATCH -e log/extval_trim-%A-%a.err
#SBATCH --array=0-95
# --array=0-95 matches the 96 samples in Mayne et al. 2020's public
# dataset (see 00_download_external_dataset.sh) - adjust to match
# sample_list.txt if you point this at a different external dataset.

module purge
module load trimgalore/0.6.10
module load fastqc/0.12.1

mapfile -t SAMPLES < "$EPICLOCK_ROOT/07_external_validation/sample_list.txt"
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID | Sample: $SAMPLE"

srun bash "$EPICLOCK_ROOT/07_external_validation/01_fastq_qc_trim.sh" "$SAMPLE"
