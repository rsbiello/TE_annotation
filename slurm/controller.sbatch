#!/usr/bin/env bash
#SBATCH --job-name=te-annotation-controller
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-controller-%j.out
#SBATCH --error=slurm-controller-%j.err

set -euo pipefail

REPOSITORY_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONFIG_FILE=${CONFIG_FILE:-config/config.yaml}
CONDA_PREFIX_DIR=${SNAKEMAKE_CONDA_PREFIX:-$REPOSITORY_DIR/.snakemake/conda}

cd "$REPOSITORY_DIR"
mkdir -p logs/slurm "$CONDA_PREFIX_DIR"

# The submission environment must already provide `snakemake` and the
# cluster-generic executor plugin. See environment.yaml.
snakemake \
    --snakefile workflow/Snakefile \
    --configfile "$CONFIG_FILE" \
    --profile profiles/slurm \
    --conda-prefix "$CONDA_PREFIX_DIR"
