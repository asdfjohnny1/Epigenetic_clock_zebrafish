#!/bin/bash
# SLURM array submission script - launches one task per sample listed in
# 07_external_validation/sample_list.txt. Adjust #SBATCH resource
# directives for your own cluster. Requires
# hpc_array_jobs/submit_01_fastq_qc_trim.sh to have completed first, and
# the danRer11 genome to already be Bismark-prepared (see
# 02_bismark_align_extract.sh header).

#SBATCH --job-name=extval_bismark
#SBATCH --partition=compute
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your-email@institution.ac.uk
#SBATCH -o log/extval_bismark-%A-%a.out
#SBATCH -e log/extval_bismark-%A-%a.err
#SBATCH --array=0-95
# --array=0-95 matches the 96 samples in Mayne et al. 2020's public
# dataset (see 00_download_external_dataset.sh) - adjust to match
# sample_list.txt if you point this at a different external dataset.

module purge
module load bismark/0.24.2
module load bowtie2/2.5.1
module load samtools/1.19

mapfile -t SAMPLES < "$EPICLOCK_ROOT/07_external_validation/sample_list.txt"
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID | Sample: $SAMPLE"

srun bash "$EPICLOCK_ROOT/07_external_validation/02_bismark_align_extract.sh" "$SAMPLE"
