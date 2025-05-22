# General
To run the algorithms, for instance `examples/compute_intersection_custom_full.jl`, do:
```bash
cd ProductPolytopesAFW
```
then set up experimental parameters in `examples/config.yml` and, finally, 
```
julia --project=.
include("examples/compute_intersection_custom_full.jl")
```

Notice that `Plots` is not installed in the `ProductPolytopesAFW` package because it is cumbersome, and we didn't want to ship our code with heavy dependencies.
For some scripts in `/examples` to run, you should have `Plots` installed in your base `Julia` environment, so that the command `using Plots` at the top of, say, `examples/compute_intersection_custom_full.jl`, won't fail

# Workflows

## Generate instances, save them as `.jld` files and run FW on them
Warning! Example scripts might be outdated
0. Decide experimental parameters in `examples/config.yml`
1. `1.1_generate_polytopes.jl` script: generate instances, find optimal solutions, then save instances to `.jld2` files
2. `1.2_compute_intersection_custom_instances.jl` script: test several FW variants on those instances

## Recommended workflow
0. Decide experimental parameters in `examples/config.yml`
1. Run either `compute_intersection_custom_full.jl` (custom instances, with parameters decided in `config.yml`) or `compute_intersection_legacy_lmos.jl` (instances defined through their `FrankWolfe.jl`'s legacy LMOs). These scripts do not save the instances to `.jld2` files

## Run on SLURM
0. Decide experimental parameters in `examples/config.yml`
1. make sure to set the 'SBATCH' parameters in the SLURM script in compliance with your SLURM server constraints
2. run the following commands
    ```
    cd /ProductPolytopesAFW
    chmod +x slurm_experiments.sh
    sbatch examples/slurm_experiments.sh examples/compute_intersection_custom_full_warmup_slurm.jl examples/results_linesearch_afw/ examples/config.yml
    ```

### Useful commands in SLURM
- `sacct -j <job_id>`: info status of one job
- `sacct -u <username>`: info status on all jobs of a user
- `sinfo -l`: general info about partitions and nodes 
- `sinfo -N -o "%N %P"`: nodes and assigned partition
- `sinfo -N -o "%N %m"`: nodes and assigned RAM
- `sinfo -N -l"`: info about nodes
- `scontrol show node=<node_name>`: shows detailed attributes about a single node
- `scancel -u <username>`: cancel all of username jobs

## Updated scripts that can be run
- `compute_inteserction_custom_full_warmup.jl`: run experiments preceded by a "warmup" experiment, so that all algorithms can be precompiled
- `slurm_experiments.sh`: run `compute_inteserction_custom_full_warmup.jl` on SLURM
- `slurm_loop.jl`: creates config files and scripts necessary to run a sequence of jobs based on `compute_inteserction_custom_full_warmup.jl` on SLURM
- `plot_from_logs.jl`: retrieves an instance from "iter_logs" and plots it
- `plot_from_logs_avg.jl`: retrieves all instance with a common pattern from "iter_logs" (for example, all instances with k=2 and n=100000) and plots average values of the trajectories in that instance