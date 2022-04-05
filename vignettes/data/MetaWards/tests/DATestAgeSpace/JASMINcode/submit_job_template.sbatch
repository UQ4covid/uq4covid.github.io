#!/bin/bash
#SBATCH -p short-serial
#SBATCH --job-name=runDesign
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --time=00:35:00
#SBATCH --array=RANGES

./runDesign.sh ${SLURM_ARRAY_TASK_ID} FILEDIR

