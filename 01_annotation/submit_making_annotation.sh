#!/bin/bash
# HPC submission script for making_annotation.R (see that file for details).
# Adjust the #SBATCH resource directives and /path/to placeholders for your own cluster.

#SBATCH --job-name=making_annotations
#SBATCH --partition=hmem
#SBATCH --qos=hmem
#SBATCH --mem=700G
#SBATCH --cpus-per-task=20
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxx
#SBATCH -o /path/to/file-%j.out
#SBATCH -e /path/to/file-%j.err


# had some issues with temporary memory so forcing a directory
# Custom temporary directory
export TMPDIR=/path/to/tmp
mkdir -p "$TMPDIR"

module purge
module load gcc/11.1.0
module add  R/4.5.1
module load OpenJDK/jdk-23.0.2

Rscript /path/to/making_annotation.R

