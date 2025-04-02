# General
To run the algorithms, for instance `examples/compute_intersection_custom_full.jl`, do:
```bash
cd BPCGProduct
```
then set up experimental parameters in `examples/config.yml` and, finally, 
```
julia --project=.
include("examples/compute_intersection_custom_full.jl")
```

Notice that `Plots` is not installed in the `BPCGProduct` package because it is cumbersome, and we didn't want to ship our code with heavy dependencies.
For some scripts in `/examples` to run, you should have `Plots` installed in your base `Julia` environment, so that the command `using Plots` at the top of, say, `examples/compute_intersection_custom_full.jl`, won't fail

# Workflows
Make sure you have Plots installed in your generic Julia environment. It is not installed upon instantiating the package, because it is very heavy and we do not want to ship it with the package

1. `1.1_generate_polytopes.jl` script: generate instances, find optimal solutions, then save instances to `.jld2` files
2. `1.2_compute_intersection_custom_instances.jl` script: test several FW variants on those instances

As an alternative:
1. `compute_intersection_custom_full.jl` (custom instances, with parameters decided in `config.yml`) or `compute_intersection_legacy_lmos.jl` (instances defined through their `FrankWolfe.jl`'s legacy LMOs). These scripts do not save the instances to `.jld2` files
