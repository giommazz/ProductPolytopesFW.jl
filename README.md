# bpcg-product
Experiments to test BPCG over the product domain on different problems

This repository contains two main components:
- `code/`: Directory containing all Julia code for experiments.
- `paper/`: Submodule linked to Overleaf for collaborative writing.

## Cloning the Repository for the First Time
To get started with the project, you need to clone the repository along with the associated Overleaf submodule. Use the following command:

```bash
git clone --recurse-submodules git@github.com:giommazz/bpcg-product.git

If you've already cloned the repository but didn't include the submodule, initialize and update the submodule with:

```bash
git submodule update --init --recursive


## Working with the Overleaf Submodule

### Pulling from, and pushing to, [OverLeaf](https://www.overleaf.com/project/65b26b4f44ada0fe2bf4be42)

To synchronize the latest changes from Overleaf into your local repository, follow these steps:

```bash
cd paper    # Navigate to `paper` (a submodule linked to OverLeaf project)
git pull overleaf master    # Pull latest changes from OverLeaf 
```

To push changes: make changes, then `add` and `commit`, and finally
```bash
git push overleaf
```

After this, changes will be visible on the OverLeaf webpage.


### Updating mother repo

Now we only have to update the mother repo

```bash
cd ..
git commit -am "update submodule reference"
git push
```
