function create_product_lmo(lmos_list, k)
    # Check if length of LMO list matches `k`
    if length(lmos_list) != k
        error("The number of LMOs provided ($(length(lmos_list))) does not match the expected number ($k).")
    end
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmos_list)
    # Create and return a ProductLMO object
    return FrankWolfe.ProductLMO(lmos_tuple)
end

# Find starting point `x0` over the product of different LMOs
function find_starting_point(prod_lmo)    
    # Prepare datafor `x0`, which is of type `FrankWolfe.BlockVector{Float64, Vector{Float64}, Tuple{Int64}}`
    # 1) Compute extreme points for each LMO in the product
    extreme_points = [FrankWolfe.compute_extreme_point(lmo, zeros(n)) for lmo in prod_lmo.lmos]    
    # 2) Generate a vector of length `k`, where each entry is a tuple `(n,)`
    block_sizes = fill((n,), length(prod_lmo.lmos))
    # 3) Compute total size of `x0`
    total_size = sum([size[1] for size in block_sizes]) 

    # Instantiate `x0`
    return FrankWolfe.BlockVector(extreme_points, block_sizes, total_size)
end