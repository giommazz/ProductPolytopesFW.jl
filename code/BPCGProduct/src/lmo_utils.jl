# `lmos.jl`
function create_product_lmo(config::Config, lmo_list)
    # Check if length of LMO list matches `k`
    if length(lmo_list) != config.k
        error("The number of LMOs provided ($(length(lmo_list))) does not match the expected number ($(config.k)).")
    end
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

function get_lmos(config::Config, filename::String; cvxhflag=false::Bool)
    # Load data and transform to Polyhedra.Polyhedron or JuMP.Model
    vertices, shifted_vertices = load_intersecting_polytopes(filename)    

    # Initialize data structures
    if cvxhflag
        lmo_list_nonintersecting = Vector{FrankWolfe.ConvexHullOracle}()    
        lmo_list_intersecting = Vector{FrankWolfe.ConvexHullOracle}()
    else
        lmo_list_nonintersecting = Vector{FrankWolfe.MathOptLMO}()
        lmo_list_intersecting = Vector{FrankWolfe.MathOptLMO}()
    end

    # Create LMOs
    for (v, vs) in zip(vertices, shifted_vertices)    
        # FrankWolfe.ConvexHullOracle objects
        if cvxhflag
            # Convert from Matrix{T} to Vector{Vector{T}}
            v = [v[i, :] for i in 1:size(v, 1)]
            vs = [vs[i, :] for i in 1:size(vs, 1)]
            # Create FrankWolfe.ConvexHullOracle objects from Matrix{T} objects ()
            lmo_nonintersecting = FrankWolfe.ConvexHullOracle(v)
            lmo_intersecting = FrankWolfe.ConvexHullOracle(vs)
            # Update data structures
            push!(lmo_list_nonintersecting, lmo_nonintersecting)
            push!(lmo_list_intersecting, lmo_intersecting)
        # FrankWolfe.MathOptLMO objects 
        else
            # Instantiate Polyedra.Polyhedron objects
            poly_nonintersecting = polytope(v)
            poly_intersecting = polytope(vs)
            # Create JuMP.Model objects from Polyedra.Polyhedron objects
            jump_poly_nonintersecting = polyhedra_to_jump(config, poly_nonintersecting)
            jump_poly_intersecting = polyhedra_to_jump(config, poly_intersecting)
            # Ensure models are optimized
            optimize!(jump_poly_nonintersecting)
            optimize!(jump_poly_intersecting)
            # Create FrankWolfe.MathOptLMO objects from JuMP.Model objects
            lmo_nonintersecting = FrankWolfe.MathOptLMO(jump_poly_nonintersecting.moi_backend)
            lmo_intersecting = FrankWolfe.MathOptLMO(jump_poly_intersecting.moi_backend)
            # Update data structures
            push!(lmo_list_nonintersecting, lmo_nonintersecting)
            push!(lmo_list_intersecting, lmo_intersecting)
        end
    end
    return lmo_list_nonintersecting, lmo_list_intersecting
end
