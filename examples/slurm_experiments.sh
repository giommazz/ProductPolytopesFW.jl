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
#SBATCH --time=14-00:00:00                      # timelimit (format is jj-hh:mm:ss). Use `sinfo` to see node time limits
#SBATCH --cpus-per-task=10                      # reserve certain amount of cores for this job
#SBATCH --mem=480G                              # total memory requested
#SBATCH -N1                                     # request a single node
#SBATCH --output=slurm-%x_%A.out                # Slurm stdout/stderr. Julia logs go to results_dir/terminal_logs
##SBATCH --partition=big                        # uncomment/adapt if your cluster requires an explicit partition
##SBATCH --nodelist=[<cluster-node-list>]       # optional: restrict to specific nodes
##SBATCH --exclude=[<nodes-to-exclude>]         # optional: exclude specific nodes
##SBATCH --mail-user=<your-email>               # optional: enable email notifications
##SBATCH --mail-type=BEGIN,END,FAIL

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

submit_dir="$(pwd)"
case "$script" in
    /*) ;;
    *) script="$submit_dir/$script" ;;
esac
case "$results_dir" in
    /*) ;;
    *) results_dir="$submit_dir/$results_dir" ;;
esac
case "$config" in
    /*) ;;
    *) config="$submit_dir/$config" ;;
esac

wrapper_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$wrapper_dir/.." && pwd)"
cd "$repo_root"

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

echo "Repository root : $repo_root"
echo "Slurm node list : ${SLURM_NODELIST:-unknown}"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "Git branch      : $(git rev-parse --abbrev-ref HEAD)"
    echo "Git commit      : $(git rev-parse HEAD)"
else
    echo "Git branch      : unavailable (not a git checkout)"
fi

# set number of threads used by Julia
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1   # force BLAS to single-threaded mode

echo "Running script: $script"
echo "Using config : $config"

config_basename=$("$JULIA_BIN" --project="$repo_root" -e '
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
srun "$JULIA_BIN" --project="$repo_root" "$script" > "$log_file" 2>&1
