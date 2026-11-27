#!/bin/bash

#SBATCH -J beta_highconf
#SBATCH -p cpu8358
#SBATCH --qos=26cores

#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00

#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

export RNA004_STAT_ROOT="${RNA004_STAT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/analysis}"



set -euo pipefail

cd ${RNA004_STAT_ROOT}/tool_vs_glori/80229_highconf_sensitivity/scripts

python3 -c "import pandas, numpy, scipy, matplotlib; print('Python packages OK')"

python3 02_highconf_stoich_sensitivity.py
