#!/bin/bash

# Run experiments on slurm servers
# How to run this:
# 1) make sure to tailor the 'SBATCH' parameters below to your directories
# 2) run the following commands
#       cd /BPCGProduct
#       chmod +x slurm_experiments.sh
#       sbatch examples/slurm_experiments.sh examples/compute_intersection_custom_full_warmup_slurm.jl examples/results_linesearch_afw/ examples/config.yml

# *************************
# Directives for the SLURM scheduler
#SBATCH --job-name=convexfeas_polytopes_afw   # job name
#SBATCH --time=14-00:00:00  # timelimit (format is jj-hh:mm:ss). Use `sinfo` to see node time limits
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --nodelist=htc-cmp502
#SBATCH --chdir=/home/htc/giommazz/bpcg-product/code/BPCGProduct/  # Navigate to dir where script you want to run is
#SBATCH --output=/home/htc/giommazz/bpcg-product/code/BPCGProduct/examples/logs/%x_%A.out # logfiles ---> %x=job name, %A=job ID
#SBATCH --partition=big  # Specify the desired partition on cluster (default: small)
##SBATCH --exclude=htc-cmp[101-148,501-532] # exclude nodes. Your job will run on nodes not in the list.

# Some useful commands to analyze partitions on SLURM
# sinfo -l: general info on cluster partitions and nodes
# sinfo -N -o "%N %P": partition that each cluster node is assigned to 
# sinfo -N -o "%N %m": memory allocated to each cluster node

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
logs_dir="examples/logs"
mkdir -p $logs_dir

echo $SLURM_NODELIST    # print node on slurm
git branch              # print git branch on which you are
git rev-parse HEAD      # print pointer to current branch

echo "Running $script with config $config"
# Extract parameters from config file name and create log file name
config_basename=$(julia --project=. -e 'using YAML, Dates; config = YAML.load_file(ARGS[1]); timestamp = Dates.format(now(), "yyyymmddHHMMSS"); oracle = config["cvxhflag"] ? "cvxho" : "lmo"; anc = config["anc_flag"] ? "anc" : "vert"; print("k", config["k"], "_n", config["n"], "_", oracle, "_", anc, "_t", timestamp)' "$config")
log_file="$logs_dir/${config_basename}_log.txt"

# Run the Julia script with the instance name and config file as parameters and redirect output to log file
srun ~/sw/julia-1.9.4/bin/julia --project=. $script > $log_file