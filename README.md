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

### Pulling Changes from Overleaf

To synchronize the latest changes from Overleaf into your local repository, follow these steps:

1. **Navigate to the `paper` directory:**
   This directory is a submodule linked to your Overleaf project.
   ```bash
   cd paper


2. Pull the latest changes from Overleaf:
```bash
git pull overleaf master  # assuming Overleaf uses 'master' as the default branch

3. Return to the root directory:
```bash
cd ..

4. Commit the updated submodule reference in your main project:
```bash
git add paper
git commit -m "Update submodule to latest Overleaf version"

5. Push the changes to your main repository:
```bash
git push origin main


### Making Local Changes and Pushing to Overleaf
If you make changes locally in the paper directory and want to push these updates back to Overleaf:

1. Navigate to the paper directory:
```bash
cd paper

2. Add and commit your changes:
```bash
git add .
git commit -m "Local updates to paper"

3. Push the changes back to Overleaf:
```bash
git push origin master  # ensure this matches Overleaf's branch name

4. Return to the root directory of your project:
```bash
cd ..

5. Commit the updated submodule reference in your main project:
```bash
git commit -am "Update submodule reference"

6. Push the changes to your main repository:
```bash
git push origin main
