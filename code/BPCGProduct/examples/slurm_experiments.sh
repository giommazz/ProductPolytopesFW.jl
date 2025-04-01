# To be executed from inside Colgen/PricingNAR: it runs `test/solomon_instance.jl` 
#   over the instances in `instances`, # with the YAML configurations in config

#!/bin/bash

# *************************
# Directives for the SLURM scheduler
#SBATCH --job-name=colgen   # job name
#SBATCH --time=07-00:00:00  # timelimit (format is jj-hh:mm:ss). Use `sinfo` to see node time limits
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --nodelist=htc-cmp019
#SBATCH --chdir=/home/htc/giommazz/ColGen-NAR/PricingNAR/  # Navigate to dir where script you want to run is
#SBATCH --output=/home/htc/giommazz/SCRATCH/log/%x_%A.out # logfiles ---> %x=job name, %A=job ID
#SBATCH --partition=big  # Specify the desired partition on cluster (default: small)
#SBATCH --exclude=htc-cmp[101-148,501-532] # exclude nodes. Your job will run on nodes not in the list.

# *************************
# Check if the correct number of arguments is passed
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 path/to/instance/dir path/to/yaml/configs/dir"
    exit 1
fi

instances_dir="$1" # Directory of instances
yaml_configs_dir="$2" # Directory of config files

# Check if the instances directory is valid
if [ ! -d "$instances_dir" ]; then
    echo "Error: The specified instances directory '$instances_dir' does not exist or is not a directory."
    exit 1
fi

# Check if the config files directory is valid
if [ ! -d "$yaml_configs_dir" ]; then
    echo "Error: The specified config files directory '$yaml_configs_dir' does not exist or is not a directory."
    exit 1
fi

# Create logs directory if it doesn't already exists (`-p` option)
logs_dir="test/logs"
mkdir -p $logs_dir
solomon_instance_script = "test/solomon_instance.jl"

echo $SLURM_NODELIST    # print node on slurm
git branch              # print git branch on which you are
git rev-parse HEAD      # print pointer to current branch

# Iterate over all XML files in the instances directory
for instance_file in "$instances_dir"/*.xml; do
    # Extract the instance name (part before .xml)
    instance_name=$(basename "$instance_file" .xml)

    # Iterate over all config files in the config files directory
    for config_file in "$yaml_configs_dir"*.yml; do
        echo "Running instance: $instance_name with config $config_file"
        
        # Extract parameters from config file name and create log file name
        config_basename=$(basename "$config_file" .yml)
        log_file="$logs_dir/${config_basename}_log.txt"

        # Run the Julia script with the instance name and config file as parameters and redirect output to log file
        srun ~/sw/julia-1.9.4/bin/julia --project=. $solomon_instance_script $instance_name $config_file > $log_file
    done
done