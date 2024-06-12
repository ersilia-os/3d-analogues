import os
import subprocess
import tempfile
import csv
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem, SDMolSupplier, SDWriter
from rdkit.Chem.rdmolfiles import MolToPDBFile

from timeout_decorator import timeout
from timeout_decorator.timeout_decorator import TimeoutError as DecoratorTimeoutError

NUM_CONF = 10

root = os.path.dirname(os.path.abspath(__file__))
tmp_dir = tempfile.mkdtemp(prefix='3danalogue-')
pdb_tmp = os.path.join(tmp_dir, "pdb_files")
sdf_tmp = os.path.join(tmp_dir, "sdf_files")

class ObabelConf():
    def __init__(self):
        pass

    @timeout(20, use_signals=False)
    def prep_single_ligand(self, smiles, id, pH, num_conformers=NUM_CONF):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        # Embed multiple conformations
        ps = AllChem.ETKDGv3()
        ps.useRandomCoords = True
        AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=ps)
        mol.SetProp('_3D', 'Yes')
        os.makedirs(pdb_tmp, exist_ok=True)
        # Output conformers to PDB
        pdb_name = os.path.join(pdb_tmp, "{}.pdb".format(id))
        MolToPDBFile(mol, pdb_name, confId=-1)
        # Protonate according to pH, convert to .sdf
        os.makedirs(sdf_tmp, exist_ok=True)
        sdf_name = os.path.join(sdf_tmp, "{}.sdf".format(id))
        cmd = "obabel {0} -pH {1} -O {2}".format(pdb_name, pH, sdf_name)
        subprocess.Popen(cmd, shell=True).wait()
        return sdf_name

    def merge_sdfs(self, out_sdfs, merged_sdf):
        mols = {}
        for s in out_sdfs:
            mol_name = s.split("/")[-1].split(".")[0]
            suppl = Chem.SDMolSupplier(s)
            for i,mol in enumerate(suppl):
                if mol is None:
                    continue
                else:
                    if mol_name not in mols:
                        mols[mol_name] = []
                    mols[mol_name].append(mol)
            #os.remove(s)
        with Chem.SDWriter(merged_sdf) as w:
            for k, v in mols.items():
                for i, mol in enumerate(v):
                    if mol is None:
                        continue
                    else:
                        mol.SetProp("_Name", "{}_conf{}".format(k,i))
                        w.write(mol)

    def generate_conformers(self, input_file, output_file):
        with open(input_file, 'r') as f:
            smiles_list = []
            id_list = []
            reader = csv.reader(f)
            headers = next(reader)
            smiles_col = None
            id_col = None
            for idx, header in enumerate(headers):
                if 'smiles' in header.lower():
                    smiles_col = idx
                elif 'id' in header.lower():
                    id_col = idx
            if smiles_col is None:
                raise ValueError("No column with 'smiles' found in the CSV file.")
            if id_col is None:
                raise ValueError("No column with 'id' found in the CSV file.")
            for row in reader:
                smiles_list.append(row[smiles_col])
                id_list.append(row[id_col])
        all_sdfs = []
        skipped_ids = []
        for smi, id in zip(smiles_list, id_list):
            try:
                print(id)
                sdf_name = self.prep_single_ligand(smi, id, 7.4)
                all_sdfs += [sdf_name]
            except DecoratorTimeoutError:
                print(f"Skipping molecule {id} due to timeout")
                skipped_ids.append(id)
                continue
        self.merge_sdfs(all_sdfs, output_file)
        shutil.rmtree(tmp_dir)