#!/bin/bash

#SBATCH -J plot_beta_strict
#SBATCH -p cpu8358
#SBATCH --qos=26cores

#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00

#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

export RNA004_STAT_ROOT="${RNA004_STAT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/analysis}"


set -euo pipefail

cd ${RNA004_STAT_ROOT}/tool_vs_glori/80229_highconf_sensitivity/scripts

python3 03_plot_strictness_figures.py
