# `lmos.jl`
function create_product_lmo(lmo_list)
    
    # Convert list of LMOs to a tuple, as required by `FrankWolfe.ProductLMO`
    lmos_tuple = Tuple(lmo_list)
    
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

# Initialize LMOs for given sets of `vertices` k=1 polytope
# Depending on `cvxhflag`, create either `FrankWolfe.ConvexHullOracle` (true) or `FrankWolfe.MathOptLMO` (false) objects.
function create_lmos(config::Config, vertices::Vector{Matrix{T}}) where T
    # Initialize data structures
    if config.cvxhflag lmo_list = Vector{FrankWolfe.ConvexHullOracle}() else lmo_list = Vector{FrankWolfe.MathOptLMO}() end

    # Create LMOs
    for v in vertices
        if config.cvxhflag # FrankWolfe.ConvexHullOracle objects
            # Convert from Matrix{T} to Vector{Vector{T}}
            v = [v[i, :] for i in 1:size(v, 1)]
            # Create FrankWolfe.ConvexHullOracle objects from Matrix{T}
            lmo = FrankWolfe.ConvexHullOracle(v)
            # Update data structures
            push!(lmo_list, lmo)
        else    # FrankWolfe.MathOptLMO objects 
            # Instantiate Polyedra.Polyhedron objects
            poly = polytope(v)
            # Create JuMP.Model objects from Polyedra.Polyhedron objects
            jump_poly = polyhedra_to_jump(config, poly)
            # Ensure models are optimized
            optimize!(jump_poly)
            # Create FrankWolfe.MathOptLMO objects from JuMP.Model objects
            lmo = FrankWolfe.MathOptLMO(jump_poly.moi_backend)
            # Update data structures
            push!(lmo_list, lmo)
        end
    end
    return lmo_list
end
# (Multiple dispatch) k polytopes
function create_lmos(config::Config, vertices::Vector{Vector{Matrix{T}}}) where T
    lmo_list = []
    for vs in vertices
        lmos = create_lmos(config, vs)
        push!(lmo_list, lmos)
    end
    return lmo_list
end
