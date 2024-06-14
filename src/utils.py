from rdkit import Chem

def extract_best_conformer(sdf_file, conformer_names, output_file):
    # Initialize a writer for the output file
    writer = Chem.SDWriter(output_file)
    # Iterate through molecules in the SDF file
    suppl = Chem.SDMolSupplier(sdf_file)
    for mol in suppl:
        if mol is None:
            continue
        mol_name = mol.GetProp("_Name")
        if mol_name in conformer_names:
            writer.write(mol)
    # Close the writer
    writer.close()