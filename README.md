# bpcg-product
Experiments to test BPCG over the product domain on different problems

This repository contains two main components:
- `code/`: containing Julia code for experiments
- `paper/`: Submodule linked to [Overleaf webpage](https://www.overleaf.com/project/65b26b4f44ada0fe2bf4be42) for collaborative writing



## Cloning the Repository for the First Time
To get started with the project, clone repository with associated OverLeaf submodule:

```bash
git clone --recurse-submodules git@github.com:giommazz/bpcg-product.git
```

If you've already cloned the repo but didn't include the submodule, initialize and update the submodule:

```bash
git submodule update --init --recursive
```



## Working with the Overleaf Submodule

### Pulling from, and pushing to, OverLeaf

Sync latest changes from remote Overleaf into local repository:

```bash
cd paper                    # `paper` is a submodule linked to OverLeaf project
git pull overleaf master    # Pull latest changes from OverLeaf 
```

Make changes, `add`, `commit`, and
```bash
git push overleaf           # changes will be visible on [OverLeaf webpage](https://www.overleaf.com/project/65b26b4f44ada0fe2bf4be42)
```

### Updating mother repo

Also push changes to remote Git repo

```bash
cd ..                       # now in main repo ../paper
git commit -am "update submodule reference"
git push
```


## Julia Setup Instructions

To set up the Julia environment necessary to run the experiments, follow these steps:

1. Clone this repository and navigate to the `code/BPCGProduct` directory.
2. Start Julia and activate the local environment
```bash
julia --project=. 
```
3. # Install and lock the exact package versions
```julia
] instantiate
```
4. Run experiments with:
```julia
include("src/BPCGProduct.jl")
BPCGProduct.runExperiment()
```
