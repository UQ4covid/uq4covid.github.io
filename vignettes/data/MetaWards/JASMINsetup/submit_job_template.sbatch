#!/bin/bash
#SBATCH -p short-serial
#SBATCH --job-name=unZIP
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --time=00:05:00
#SBATCH --array=RANGES

./unZIP.sh ${SLURM_ARRAY_TASK_ID} FILEDIR

