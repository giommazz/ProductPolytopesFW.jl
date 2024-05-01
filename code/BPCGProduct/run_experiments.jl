# Assuming your module's code is in 'src/BPCGProduct.jl'
include("src/BPCGProduct.jl")  # Include the module
using .BPCGProduct  # Use the module

# Call the main function to run your experiments
BPCGProduct.main()

