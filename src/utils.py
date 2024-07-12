import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

class SDFProcessor():
    def __init__(self,output_folder):
        self.output_folder = output_folder

    def extract_best_conformer(self, input_file,  output_file, conformer_names,):
        writer = Chem.SDWriter(os.path.join(self.output_folder,output_file))
        suppl = Chem.SDMolSupplier(os.path.join(self.output_folder,input_file), removeHs=False)
        for i,mol in enumerate(suppl):
            if mol is not None:
                mol_name = mol.GetProp("_Name")
                if mol_name in conformer_names:
                    writer.write(mol)
            else:
                print(f"Warning: Molecule at index {i} is invalid and will be skipped.")
        # Close the writer
        writer.close()

    def get_original_names(self, input_file):
        original_names = []
        suppl = Chem.SDMolSupplier(os.path.join(self.output_folder, input_file), removeHs=False)
        for i,mol in enumerate(suppl):
            if mol is not None:
                mol_id = mol.GetProp("_Name")
                original_names += [mol_id]
            else:
                print(f"Warning: Molecule at index {i} is invalid and will be skipped.")
        return original_names

    def append_names_to_sdf(self, input_file, output_file, original_names):
        suppl = Chem.SDMolSupplier(os.path.join(self.output_folder, input_file),removeHs=False)
        writer = Chem.SDWriter(os.path.join(self.output_folder,output_file))
        for i, (mol, mol_id) in enumerate(zip(suppl, original_names)):
            if mol is not None:
                mol.SetProp("_Name", mol_id)
                writer.write(mol)
            else:
                print(f"Warning: Molecule at index {i} is invalid and will be skipped.")
        writer.close()

    def add_pose_to_name(self, input_file, output_file):
        writer = Chem.SDWriter(os.path.join(self.output_folder, output_file))
        suppl = Chem.SDMolSupplier(os.path.join(self.output_folder, input_file), removeHs=False)
        conf2pose = {}
        for mol in suppl:
            if mol is None:
                continue
            mol_name = mol.GetProp("_Name")
            if mol_name in conf2pose.keys():
                i_ = len(conf2pose[mol_name])
                conf2pose[mol_name].append(mol_name+"_"+f"pose{i_}")
                mol.SetProp("_Name", mol_name+"_"+f"pose{i_}")
                writer.write(mol)
            else:
                conf2pose[mol_name] = [mol_name+"_"+"pose0"]
                mol.SetProp("_Name", mol_name+"_"+"pose0")
                writer.write(mol)
        writer.close()

    def add_hs(self, input_file, output_file):
        supplier = Chem.SDMolSupplier(os.path.join(self.output_folder, input_file), removeHs=False)
        mols = [mol for mol in supplier if mol is not None]
        mols_with_hs = []
        for mol in mols:
            # Add hydrogens without moving other atoms considerably
            mol_with_hs = Chem.AddHs(mol, addCoords=True)
            # Get the atom indices of the added hydrogens
            hydrogen_indices = [atom.GetIdx() for atom in mol_with_hs.GetAtoms() if atom.GetAtomicNum() == 1]
            # Create a force field for the molecule with hydrogens
            ff = AllChem.UFFGetMoleculeForceField(mol_with_hs, confId=0)
            # Constrain the positions of non-hydrogen atoms
            for i in range(mol_with_hs.GetNumAtoms()):
                if i not in hydrogen_indices:
                    ff.AddFixedPoint(i)
            # Optimize the molecule, allowing only hydrogen atoms to move
            ff.Minimize(maxIts=200)
            mols_with_hs += [mol_with_hs]
        w = Chem.SDWriter(os.path.join(self.output_folder, output_file))
        for mol in mols_with_hs:
            w.write(mol)
        w.close()

    def create_matrix_csv(self, original_names):
        with open(os.path.join(self.output_folder,'matrix-sensaas.txt'), 'r') as file:
            line = file.readline().strip()
        values = line.split()
        df = pd.DataFrame([values])
        df = df.transpose()
        df = df.rename(columns={0:"sensaas_score"})
        df["conf_id"] = original_names
        df["id"] = [original_name.split("_")[0] for original_name in original_names]
        df["idx_pose"] = df.index #save the original index position
        df["idx_pose"] = df.apply(lambda row: (row["idx_pose"]+1), axis=1)
        df.to_csv(os.path.join(self.output_folder, "matrix_sensaas.csv"), index=False)
        df['sensaas_score'] = df['sensaas_score'].astype(float)
        best_df = df.loc[df.groupby('id')['sensaas_score'].idxmax()] #this reorders the molecules, so the .csv will not be in the same order as the sdf, use the idx column
        best_df.to_csv(os.path.join(self.output_folder, "matrix_sensaas_best.csv"), index=False)
        return best_df