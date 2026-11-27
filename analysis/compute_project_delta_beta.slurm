#!/usr/bin/env bash
#SBATCH --job-name=delta_beta
#SBATCH -p cpu8358
#SBATCH --qos=26cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --output=delta_beta.%j.out
#SBATCH --error=delta_beta.%j.err

export RNA004_STAT_ROOT="${RNA004_STAT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/analysis}"


set -euo pipefail

SCRIPT="${RNA004_STAT_ROOT}/compute_project_delta_beta.py"

echo "Job started: $(date)"
echo "Host: $(hostname)"

python3 "$SCRIPT" \
    --project all \
    --chunk-size 500000

echo "Job completed: $(date)"`
