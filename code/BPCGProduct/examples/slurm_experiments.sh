# Run experiments on slurm servers
# How to run this:
# 1) make sure to tailor the 'SBATCH' parameters below to your directories
# 2) run the following commands
#       cd /BPCGProduct
#       chmod +x slurm_experiments.sh
#       ./examples/slurm_experiments.sh examples/compute_intersection_warmup.jl examples/compute_intersection_custom_full.jl examples/results_linesearch_afw/ examples/config.yml



#!/bin/bash

# *************************
# Directives for the SLURM scheduler
#SBATCH --job-name=convexfeas_polytopes_afw   # job name
#SBATCH --time=07-00:00:00  # timelimit (format is jj-hh:mm:ss). Use `sinfo` to see node time limits
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --nodelist=htc-cmp019
#SBATCH --chdir=/home/htc/giommazz/bpcg-product/code/BPCGProduct/examples  # Navigate to dir where script you want to run is
#SBATCH --output=/home/htc/giommazz/SCRATCH/log/%x_%A.out # logfiles ---> %x=job name, %A=job ID
#SBATCH --partition=big  # Specify the desired partition on cluster (default: small)
#SBATCH --exclude=htc-cmp[101-148,501-532] # exclude nodes. Your job will run on nodes not in the list.

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
config_basename=$(basename "$config" .yml)
log_file="$logs_dir/${config_basename}_log.txt"

# Run the Julia script with the instance name and config file as parameters and redirect output to log file
srun ~/sw/julia-1.9.4/bin/julia --project=. $script > $log_file