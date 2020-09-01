#!/bin/bash
#SBATCH -p short-serial
#SBATCH --job-name=createSummary
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=RANGES
#SBATCH --time 02:00:00
#SBATCH --time-min 00:20:00

module add jaspy

./createSum.sh ${SLURM_ARRAY_TASK_ID}

