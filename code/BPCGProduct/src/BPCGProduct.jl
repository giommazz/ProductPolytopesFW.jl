module BPCGProduct

using FrankWolfe

# Example function that utilizes FrankWolfe
function runExperiment()
    
    f(p) = sum(abs2, p)  # objective function
    grad!(storage, p) = storage .= 2p  # in-place gradient computation
    # function d ⟻ argmin ⟨p,d⟩ st. p ∈ Δ
    lmo = FrankWolfe.ProbabilitySimplexOracle(1.)
    p0 = [1., 0., 0.]
    
    println("Running experiment with Frank Wolfe...")
    
    p_opt, _ = frank_wolfe(f, grad!, lmo, p0; verbose=true);
end


end # module