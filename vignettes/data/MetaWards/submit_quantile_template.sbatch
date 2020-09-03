#!/bin/bash
#SBATCH -p par-single
#SBATCH --job-name=createSummaryR
#SBATCH -n 10
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=RANGES

module add jaspy

./createQuantiles.sh ${SLURM_ARRAY_TASK_ID} 20 10

