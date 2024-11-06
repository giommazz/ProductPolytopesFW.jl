# Workflow

1. run `1.1_generate_polytopes.jl`: generate instances, find optimal solutions, then save instances to `.jld2` files
2. run `1.2_compute_intersection_custom_instances.jl`: test several FW variants on those instances

As an alternative, the scripts `compute_intersection_custom_full.jl`, or `compute_intersection_legacy_lmos.jl` do the two steps above (with custom instances, or using the FW legacy LMOs) but don't save the instances to `.jld2` files

Note: to decide stepsizes, go into `product_algorithms.jl`