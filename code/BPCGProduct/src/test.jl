module test
using FrankWolfe

k = 10
n = 10
target_tolerance = 1e-9
max_iterations = 10000

# Function to compute the pairwise distance objective
function f(x)
    sum_dist = 0.0
    for i in 1:k    
        for j in i+1:k
            xi = x.blocks[i]
            xj = x.blocks[j]
            curr = sum((xi - xj).^2)
            sum_dist += curr
        end
    end
    return 0.5 * sum_dist
end

# TODO: SEE BLOCKVECTOR STRUCTURE AND PROD_LMO STRUCTURE

# Gradient computation for tuple of vectors
function grad!(storage, x)
    #println("°°°°°°°° storage (of type $(typeof(storage))): ", storage, typeof(storage))
    #println("°°°°°°°°x (of type $(typeof(x))): ", x)
    for i = 1:k
        #println("°°°°°°°° xᵢ (of type $(typeof(x.blocks[i]))): ", x.blocks[i])
        sum_terms = zeros(n)
        for j = 1:k
            #println("°°°°°°°° xᵢ (of type $(typeof(x.blocks[j]))): ", x.blocks[j])
            if i != j
                sum_terms .+= x.blocks[j]
            end
        end
        storage.blocks[i] .= 0.5 * ((k-1) * x.blocks[i] - sum_terms)
    end
end


function create_block_vector(lmos)
    # Define the ProductLMO
    prod_lmo = FrankWolfe.ProductLMO(lmos)
    
    # Compute extreme points for each LMO in the ProductLMO
    extreme_points = [FrankWolfe.compute_extreme_point(lmo, zeros(n)) for lmo in lmos]
    
    # Prepare sizes for the BlockVector
    block_sizes = fill((n,), length(lmos))
    total_size = sum([size[1] for size in block_sizes])  # Calculate total size

    # Create the BlockVector (FrankWolfe.BlockVector{Float64, Vector{Float64}, Tuple{Int64}})
    x0 = FrankWolfe.BlockVector(extreme_points, block_sizes, total_size)
    
    return x0
end

# Run ALM
function run_FW(order, prod_lmo, x0, k, n, target_tolerance, max_iterations)   
    trajectories = []

    # Check (https://zib-iol.github.io/FrankWolfe.jl/dev/examples/docs_10_alternating_methods/): performing a full (simulatenous) BPCG update at each iteration, 
    # 	by running `alternating_linear_minimization` with `blended_pairwise_conditional_gradient` inside
    x, v, primal, dual_gap, trajectory_data = FrankWolfe.block_coordinate_frank_wolfe(
        f,
        grad!,
        prod_lmo,
        x0,
        update_order=order,
        epsilon=target_tolerance,
        max_iteration=max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
# Run BPCG
function run_FW(prod_lmo, x0, k, n, target_tolerance, max_iterations)   
    trajectories = []
    # Check (https://zib-iol.github.io/FrankWolfe.jl/dev/examples/docs_10_alternating_methods/): performing a full (simulatenous) BPCG update at each iteration, 
    # 	by running `alternating_linear_minimization` with `blended_pairwise_conditional_gradient` inside
    x, v, primal, dual_gap, _, trajectory_data = FrankWolfe.blended_pairwise_conditional_gradient(
        f,
        grad!,
        prod_lmo,
        x0,
        epsilon=target_tolerance,
        max_iteration=max_iterations,
        line_search=FrankWolfe.Shortstep(one(Int)),
        print_iter=max_iterations / 10,
        memory_mode=FrankWolfe.InplaceEmphasis(),
        verbose=true,
        trajectory=true,
    );
    push!(trajectories, trajectory_data)    
    return trajectories
end
	
# Function to retrieve an LMO by name using eval
function get_lmo(lmo_name)
    lmo = eval(Symbol(lmo_name))
    return lmo
end

function create_lmos(lmo_names::Vector{String})
    lmos = Tuple(get_lmo(name) for name in lmo_names)
    return lmos
end

# Function to generate LMO names based on the number k
function generate_lmo_names(k)
    lmo_names = ["lmo$i" for i in 1:k]
    return lmo_names
end


# Main execution function
function main()

	# Setup Linear Minimization Oracles for the polytopes
    lmo1 = FrankWolfe.ProbabilitySimplexOracle(1.0)
    lmo2 = FrankWolfe.ProbabilitySimplexOracle(123.0)
    lmo3 = FrankWolfe.ProbabilitySimplexOracle(123.0)
    lmo4 = FrankWolfe.ProbabilitySimplexOracle(123.0)
    lmo5 = FrankWolfe.ProbabilitySimplexOracle(123.0)
    lmo6 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))
    lmo7 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))
    lmo8 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))
    lmo9 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))
    lmo10 = FrankWolfe.ScaledBoundLInfNormBall(-ones(n), ones(n))

    
    # TODO: VERIFICA CHE LE TRE RIGHE QUI SOTTO FUNZIONINO
    #lmo_names = generate_lmo_names(k)
    #lmos = create_lmos(lmo_names)
    lmos = (lmo1, lmo2, lmo3, lmo4, lmo5, lmo6, lmo7, lmo8, lmo9, lmo10)

    prod_lmo = FrankWolfe.ProductLMO(lmos)

    # TODO: CHECK src/alternating_linear_minimization FrankWolfe.jl
    #x0 = FrankWolfe.BlockVector([-ones(n), [i == 1 ? 1 : 0 for i in 1:n]], fill((n,), k), k * n)
    #x1 = compute_extreme_point(FrankWolfe.ProductLMO(lmos), tuple(fill(zeros(n), k)...))
    x0 = create_block_vector(lmos)
    #println("°°°°°°°° x₀ (type $(typeof(x0))): ", x0)
    #println("°°°°°°°° x₁ (type $(typeof(x1))): ", x1)

    #x0 = tuple([FrankWolfe.compute_extreme_point(lmo, zeros(n)) for lmo in lmos]...)

    # Run BPCG: FrankWolfe.blended_pairwise_conditional_gradient, 
    println("\n\n---------------------------------------------------------")
    println("Blended Pairwise Conditional Gradient")
    bpcg_trajectories = run_FW(prod_lmo, x0, k, n, target_tolerance, max_iterations)

    # Run ALM with block-coordinate Frank-Wolfe: FrankWolfe.block_coordinate_frank_wolfe
    println("\n\n---------------------------------------------------------")
    println("Alternating Linear Minimization")
    alm_trajectories = run_FW(FrankWolfe.CyclicUpdate(), prod_lmo, x0, k, n, target_tolerance, max_iterations)

end

end # module BPCGProduct