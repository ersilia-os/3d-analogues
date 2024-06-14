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
We provide an example using data from Mediouni et al, 2019. Using DCA as a query molecule, and a list of 12 molecules (DCA and 11 analogues), we try the docking to the Tat1 HIV protein (1k5k_1.pdbqt from PDBank). 

As a result, you will obtain the docking scores and 3D similarity to DCA.

```bash
cd 3d-analogues
python src/main.py -q example/dca.csv -s example/mediouni2019.csv -o results -cdpkit True
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
We use the VSFlow pipeline (please cite [Jung et al, 2023](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00703-1)if you use it) To calculate the following metrics of overlap between a query molecule and the list of smiles. If no query molecule is provided, the scorer will not run. To facilitate the analysis, we only keep the best conformer according to the docking for screening. If you already have a preferred conformer for the query molecule, you can input that directly as an sdf file and it will be used. 
* ComboScore: average of Shape Similarity and 3D Fingerprint Similarity
* Shape Similarity: shape-based similarity (using RDKIT Open3DAlign)
* FP3D Similarity: similarity calculated with 3D pharmacopore fingerprints (RDKIT Pharm2D)