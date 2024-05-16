# `lmos.jl`
function create_product_lmo(config::Config, lmos_list)
    
    println(typeof(lmos_list))
    println("debug in lmos_list.jl and remember the readline()!!!")
    readline()
    
    # Check if length of LMO list matches `k`
    if length(lmos_list) != config.k
        error("The number of LMOs provided ($(length(lmos_list))) does not match the expected number ($(config.k)).")
    end
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmos_list)
    # Create and return a ProductLMO object
    return FrankWolfe.ProductLMO(lmos_tuple)
end

# Find starting point `x0` over the product of different LMOs
function find_starting_point(config::Config, prod_lmo::FrankWolfe.ProductLMO)
    # Prepare datafor `x0`, which is of type `FrankWolfe.BlockVector{Float64, Vector{Float64}, Tuple{Int64}}`
    # 1) Compute extreme points for each LMO in the product
    extreme_points = [FrankWolfe.compute_extreme_point(lmo, zeros(config.n)) for lmo in prod_lmo.lmos]
    # 2) Generate a vector of length `k`, where each entry is a tuple `(n,)`
    block_sizes = fill((config.n,), length(prod_lmo.lmos))
    # 3) Compute total size of `x0`
    total_size = sum([size[1] for size in block_sizes])

    # Instantiate `x0`
    return FrankWolfe.BlockVector(extreme_points, block_sizes, total_size)
end