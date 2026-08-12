#!/bin/bash
# SLURM array submission script - launches one task per scaffold file
# listed in scaffold_file_list.txt (produced by the step 1-5 script).
# Adjust #SBATCH resource directives for your own cluster.

#SBATCH --job-name=EN_S6_cluster_array
#SBATCH --partition=hmem
#SBATCH --qos=hmem
#SBATCH --mem=300G
#SBATCH --cpus-per-task=3
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your-email@institution.ac.uk
#SBATCH -o $PROJECT_ROOT/multitissue/making_clocks/R/log/cluster/EN_S6_cluster_array-%A-%a.out
#SBATCH -e $PROJECT_ROOT/multitissue/making_clocks/R/log/cluster/EN_S6_cluster_array-%A-%a.err


module purge
module load gcc/11.1.0
module load  R/4.5.1
module load OpenJDK/jdk-23.0.2


# Read the scaffold file list into an array
mapfile -t FILE_LIST < $PROJECT_ROOT/multitissue/making_clocks/R/results/cluster_clock_new/step3.2/scaffold_file_list.txt

# Get the file corresponding to this array task
INPUT_FILE=${FILE_LIST[$SLURM_ARRAY_TASK_ID]}

echo "SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID"
echo "Processing file: $INPUT_FILE"


srun Rscript $PROJECT_ROOT/multitissue/making_clocks/R/scripts/array/post_evolution/new/Cluster_step6_elastic_net_final_new.R $INPUT_FILE

