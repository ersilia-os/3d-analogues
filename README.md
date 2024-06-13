# 3D Analogues
This package provides several evaluation metrics for identifying the best analogues of molecules in the 3D space.

# Installation

Create a conda environment and install the packages listed in requirements.txt using pip:
```bash
conda create -n analogues python=3.11
conda activate analogues
pip install -r requirements.txt
```

Install the open source version of PyMol:
```bash
conda install -c conda-forge pymol-open-source
```

If the conformer generation is done with OpenBabel, the user wil need to install OpenBabel in the system.

# Usage

```bash
cd 3d-analogues
python src/main.py -q example/query_mol.csv -s example/molecules.csv -o results -cdpkit True
```

# How it works
Using a starting molecule and a list of putative analogue candidates (we recommend looking at [ChemSampler](https://github.com/chem-sampler)) it will provide a numeric score based on twe metrics: docking to protein of interest and 3D colocalisation with query molecule

## 1. Conformer generation
The first step is to convert the SMILES of the molecules into 3D conformers. We can do so with:
* OpenBabel
* CDPKit (flag -cdpkit True) - uses code from this [repository](https://github.com/ersilia-os/smiles-to-3d)
In both cases it defaults to generating 10 conformers per molecule with the minimum energy.

## 2. Docking
This package attempts the docking to a protein of interest. The following files are required in the `proteins` folder:
* `protein.pdbqt`: protein file from pdbank, for example. The example is the 1k5k_1 TAT Protein (HIV)
* `box_coords.json`: JSON file containing the x,y,z coordinates of the binding box and its size. An example is provided.

Optional: if a file named `residue_coords.json` is found in the `proteins` folder, the distance between the molecule and the selected protein residues will be calculated and used to decide the best conformer and pose of each docked molecule (`best_docking_results.csv`). Currently the docking score and the distance score weight 50% each. The results for all conformers (up to 10) and poses (up to 10) are stored in the indicated results folder under `all_docking_results.csv`.

## 3. 3D shape scorer
We use the VSFlow pipeline (please cite [Jung et al, 2023](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00703-1)if you use it) To calculate the following metrics of overlap between a query molecule and the list of smiles. If no query molecule is provided, the scorer will not run.
* ComboScore
* Shape Similarity
* FP3D Similarity 
