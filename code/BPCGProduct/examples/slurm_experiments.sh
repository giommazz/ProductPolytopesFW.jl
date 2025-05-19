#!/bin/bash
# Run experiments on slurm servers
# How to run this:
# 1) make sure to tailor the 'SBATCH' parameters below to your directories
# 2) run the following commands
#       cd /ProductPolytopesAFW
#       chmod +x slurm_experiments.sh
#       sbatch examples/slurm_experiments.sh examples/compute_intersection_custom_full_warmup_slurm.jl examples/results_linesearch_afw/ examples/config.yml

# *************************
# Directives for the SLURM scheduler
#SBATCH --job-name=convexfeas_polytopes_afw   # job name
#SBATCH --time=14-00:00:00      # timelimit (format is jj-hh:mm:ss). Use `sinfo` to see node time limits
#SBATCH --cpus-per-task=2       # reserve a certain amount of cores for this job
#SBATCH --mem=480G
#SBATCH -N1 --nodelist=htc-cmp[501-532] # choose one node from the list in the square brackets
#SBATCH --chdir=/home/htc/giommazz/afw-product/code/ProductPolytopesAFW/  # Navigate to dir where script you want to run is
#SBATCH --output=/home/htc/giommazz/afw-product/code/ProductPolytopesAFW/examples/logs/%x_%A.out # logfiles ---> %x=job name, %A=job ID
#SBATCH --partition=big  # Specify the desired partition on cluster (default: small)
##SBATCH --exclude=htc-cmp[101-148,501-532] # exclude nodes. Your job will run on nodes not in the list.

# ensures the script halts on any error, so you don't accidentally skip the trap below
set -euo pipefail

# Print start time
echo "Job started at: $(date)"

# Ensure we always log an end time (unless the shell is forcibly killed): 
#       this registers the `echo` command to run when the bash script exits, for any exit path other than untrappable signals 
#       Even if `srun` returns a nonzero status because the job was OOM-killed, bash still invokes the EXIT trap, printing your end time.
trap 'echo "Job ended at:   $(date)"' EXIT ERR

# *************************
# Check if the correct number of arguments is passed: the user must specify a directory to save the results
if [ "$#" -lt 2 ]; then
    echo "Usage: $0     path/to/script     path/to/results/dir     path/to/config/file"
    exit 1
fi

script="$1" # Script to be run
results_dir="$2" # Directory containing the results
config="$3" # YAML config file containing params for running experiments

# Preeliminary controls
if [ ! -f "$script" ]; then
    echo "Error: The specified script '$script' does not exist."
    exit 1
fi
if [ ! -d "$results_dir" ]; then
    echo "Error: The specified path '$results_dir' does not exist or is not a directory."
    exit 1
fi
if [ ! -f "$config" ]; then
    echo "Error: The specified file '$config' does not exist."
    exit 1
fi
if [[ "$config" != *.yaml && "$config" != *.yml ]]; then
    echo "Error: The specified file '$config' is not a YAML file."
    exit 1
fi

# Create logs directory if it doesn't already exists (`-p` option)
logs_dir="examples/results_linesearch_afw/logs"
mkdir -p $logs_dir

echo $SLURM_NODELIST    # print node on slurm
git branch              # print git branch on which you are
git rev-parse HEAD      # print pointer to current branch

# set number of threads used by Julia
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=1   # force BLAS to single‑threaded mode

echo "Running $script with config $config"
# Extract parameters from config file name and create log file name
config_basename=$(julia --project=. -e '
    using YAML, Dates;
    config = YAML.load_file(ARGS[1]);
    timestamp  = Dates.format(now(), "yyyymmdd_HHMMSS");
    oracle = config["cvxhflag"] ? "cvxho" : "lmo";
    anc    = config["anc_flag"] ? "anc"   : "vert";
    println("k",  config["k"],
            "_n", config["n"],
            "_i", config["max_iterations"],
            "_s", config["seed"],
            "_",  oracle,
            "_",  anc,
            "_t", timestamp)
' "$config")

log_file="$logs_dir/${config_basename}_log.txt"

# Run the Julia script with the instance name and config file as parameters and redirect output to log file
srun ~/sw/julia-1.9.4/bin/julia --project=. $script > $log_file