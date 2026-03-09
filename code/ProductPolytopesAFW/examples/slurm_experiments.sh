#!/bin/bash
# Run one experiment script on SLURM with a specific config file.
#
# Usage:
#   sbatch examples/slurm_experiments.sh \
#     examples/compute_intersection_point_clouds_full_warmup.jl \
#     examples/results_linesearch_point_clouds \
#     examples/config.yml

# *************************
# Directives for SLURM scheduler
#SBATCH --job-name=prod_polytope_intersect_fw   # job name
#SBATCH --time=14-00:00:00      # timelimit (format is jj-hh:mm:ss). Use `sinfo` to see node time limits
#SBATCH --cpus-per-task=10      # reserve a certain amount of cores for this job
#SBATCH --mem=480G
#SBATCH -N1 --nodelist=htc-cmp[501-532] # choose one node from the list in the square brackets
#SBATCH --chdir=/home/htc/giommazz/product-polytopes-afw/code/ProductPolytopesAFW/  # Navigate to dir where script you want to run is
#SBATCH --output=/home/htc/giommazz/product-polytopes-afw/code/ProductPolytopesAFW/examples/logs/%x_%A.out # logfiles ---> %x=job name, %A=job ID
#SBATCH --partition=big  # Specify the desired partition on cluster (default: small)
##SBATCH --exclude=htc-cmp[101-148,501-532] # exclude nodes. Your job will run on nodes not in the list.

set -euo pipefail

echo "Job started at: $(date)"
trap 'echo "Job ended at:   $(date)"' EXIT ERR

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <path/to/script.jl> <path/to/results_dir> <path/to/config.yml>"
    exit 1
fi

script="$1"
results_dir="$2"
config="$3"
JULIA_BIN="${JULIA_BIN:-julia}"

# Preliminary checks
if [ ! -f "$script" ]; then
    echo "Error: The script '$script' does not exist."
    exit 1
fi
if [ ! -f "$config" ]; then
    echo "Error: The config file '$config' does not exist."
    exit 1
fi
if [[ "$config" != *.yaml && "$config" != *.yml ]]; then
    echo "Error: '$config' is not a YAML file."
    exit 1
fi

mkdir -p "$results_dir"
logs_dir="$results_dir/terminal_logs"
mkdir -p "$logs_dir"

echo "$SLURM_NODELIST"   # print node on slurm
git branch               # print git branch on which you are
git rev-parse HEAD       # print pointer to current branch

# set number of threads used by Julia
export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS=1   # force BLAS to single-threaded mode

echo "Running script: $script"
echo "Using config : $config"

config_basename=$("$JULIA_BIN" --project=. -e '
    using YAML, Dates
    config = YAML.load_file(ARGS[1])
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    oracle = get(config, "cvxhflag", false) ? "cvxh" : "mathopt"
    backend = string(get(config, "convex_hull_backend", "na"))
    intersection_cfg = get(config, "intersection", Dict{Any,Any}())
    anchor = string(get(intersection_cfg, "anchor", "na"))
    sanitize(s) = replace(lowercase(String(s)), r"[^a-z0-9_-]+" => "-")
    println(
        "k", config["k"],
        "_n", config["n"],
        "_i", config["max_iterations"],
        "_s", config["seed"],
        "_", oracle,
        "_", sanitize(backend),
        "_", sanitize(anchor),
        "_t", timestamp
    )
' "$config")

script_basename="$(basename "$script" .jl)"
log_file="$logs_dir/${script_basename}_${config_basename}.log"

echo "Julia executable: $JULIA_BIN"
echo "Terminal log    : $log_file"

# Run Julia script and save stdout/stderr to terminal log
srun "$JULIA_BIN" --project=. "$script" > "$log_file" 2>&1