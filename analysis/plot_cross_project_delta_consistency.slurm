#!/usr/bin/env bash
#SBATCH --job-name=cross_proj_delta
#SBATCH -p cpu8358
#SBATCH --qos=26cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=1-00:00:00
#SBATCH --output=cross_proj_delta.%j.out
#SBATCH --error=cross_proj_delta.%j.err

export RNA004_STAT_ROOT="${RNA004_STAT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/analysis}"


set -euo pipefail

SCRIPT="${RNA004_STAT_ROOT}/plot_cross_project_delta_consistency.py"

echo "Job started: $(date)"
echo "Host: $(hostname)"

python3 "$SCRIPT"

echo "Job completed: $(date)"
