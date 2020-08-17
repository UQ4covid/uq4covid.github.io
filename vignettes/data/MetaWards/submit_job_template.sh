#!/bin/bash
#BSUB -q par-multi
#BSUB -J createSummary
#BSUB -n 20
#BSUB -o %J.o
#BSUB -e %J.e
#BSUB -W 10:00

module load jaspy

run-multi <<EOF
REPLACE_STUFF
EOF

