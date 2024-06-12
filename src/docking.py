import os
import subprocess
import sys
import pandas as pd
from scipy.stats import rankdata
import json
from rdkit import Chem
from rdkit.Chem import PandasTools


root = os.path.dirname(os.path.abspath(__file__))
SMINA = os.path.join(root, "..","smina", "smina.static") 
RECEPTOR = os.path.join(root, "..", "proteins", "protein.pdbqt")
RESIDUE_COORDS =  os.path.join(root, "..", "proteins", "residue_coords.json")
BOX_COORDS =  os.path.join(root, "..", "proteins", "box_coords.json")


class ProteinDocker():
    def __init__(self, output_folder):
        self.output_folder = output_folder
        self.docking_result = os.path.join(self.output_folder, "docking_output.sdf")


    def _read_box_coords(self):
        if os.path.exists(BOX_COORDS):
            with open(BOX_COORDS, "r") as json_file:
                data = json.load(json_file)
                if "XYZ" in data and "BOX" in data:
                    return data["XYZ"], data["BOX"]
                else:
                    raise ValueError("JSON file does not contain 'XYZ' or 'BOX' parameters.")
        else:
            raise FileNotFoundError("box_coords.json file not found, and is required")

    def _docking_cmd(self):
        xyz, box = self._read_box_coords()
        center_x = xyz[0]
        center_y = xyz[1]
        center_z = xyz[2]
        ligands = os.path.join(self.output_folder, "conformers.sdf")
        log = os.path.join(self.output_folder, "log.txt")
        cmd = (
            SMINA
            + " -r "
            + RECEPTOR
            + " -l "
            + ligands
            + " -o "
            + self.docking_result
            + " --log "
            + log
            + f" --center_x {center_x}"
            + f" --center_y {center_y}"
            + f" --center_z {center_z}"
            + f" --size_x {box} "
            + f" --size_y {box} "
            + f" --size_z {box} "
            + " --exhaustiveness 10 "
            + " --num_modes 9 "
            + " --addH off "
            )
        subprocess.Popen(cmd, shell=True).wait()
    
    def _load_docking_results(self):
        df = PandasTools.LoadSDF(self.docking_result, embedProps=True, molColName=None)
        return df
    
    def get_docking_score(self):
        self._docking_cmd()
        df = self._load_docking_results()
        pose_counter = {}
        new_ids = []
        for _, row in df.iterrows():
            id = row['ID']
            if id not in pose_counter:
                pose_counter[id] = 0
            new_id = f"{id}_pose{pose_counter[id]}"
            pose_counter[id] += 1
            new_ids.append(new_id)
        df['ID'] = new_ids
        df["id"] = [new_id.split("_")[0] for new_id in new_ids]
        df.rename(columns={"minimizedAffinity": "docking_score", "ID": "conf_id"}, inplace=True)
        df["docking_rank"] = len(df) + 1 - rankdata(df["docking_score"]) # we want values more negative
        df = df[["id", "conf_id", "docking_score", "docking_rank"]]
        return df

class ResidueDistanceCalc(ProteinDocker):
    def __init__(self, output_folder):
        super().__init__(output_folder)

    def _get_residue_coords(self):
        with open(RESIDUE_COORDS, 'r') as json_file:
            residue_coords_dict = json.load(json_file)
        return residue_coords_dict

    def _extract_atom_coords_from_sdf(self, sdf_file):
        sdf_supplier = Chem.SDMolSupplier(sdf_file)
        positions = {}
        for mol in sdf_supplier:
            if mol is None:
                continue
            if not mol.HasProp("_Name"):
                raise ValueError("Molecule does not have a name (_Name property).")
            # Store the coordinates in the positions dictionary
            mol_title = mol.GetProp("_Name")
            coords = []
            # Iterate through each atom in the molecule
            for atom in mol.GetAtoms():
                # Get the atom index
                atom_idx = atom.GetIdx()
                # Get the 3D coordinates of the atom
                pos = mol.GetConformer().GetAtomPosition(atom_idx)
                x, y, z = pos.x, pos.y, pos.z
                atom_type = atom.GetSymbol()
                coords.append((atom_type, (x, y, z)))
            if mol_title in positions:
                positions[mol_title] += [{i: (atom_type, (x, y, z)) for i, (atom_type, (x, y, z)) in enumerate(coords, 1)}]
            else:
                positions[mol_title] = [{i: (atom_type, (x, y, z)) for i, (atom_type, (x, y, z)) in enumerate(coords, 1)}]
        return positions

    def _calculate_distance(self, selected_coords, target_coords):
        # Unpack coordinates for the selected residue
        x1, y1, z1 = selected_coords
        # Unpack coordinates for the target residue
        x2, y2, z2 = target_coords
        # Calculate Euclidean distance
        distance = ((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)**0.5
        return distance
    
    def _keep_heteroatoms_only(self):
        positions = self._extract_atom_coords_from_sdf(self.docking_result)
        positions_ha = {}
        for k,v in positions.items():
            ha_list = []
            for i in range(len(v)):
                heteroatoms = {k_: v_ for k_,v_ in v[i].items() if v_[0] not in ["C1", "C2", "C3", "C4"]}
                ha_list += [heteroatoms]
            positions_ha[k] = ha_list
        return positions_ha
    
    def _get_distance_by_pose(self):
        positions_ha = self._keep_heteroatoms_only()
        residue_coords_dict = self._get_residue_coords()
        #calculate distances to relevant atoms in residues
        all_distances = {}
        for k,v in positions_ha.items():
            dist_by_pose = {}
            for i in range(len(v)):
                dist_by_res = {}
                for res, res_coord in residue_coords_dict.items():
                    distances = []
                    for k_, v_ in v[i].items():
                        dist = self._calculate_distance(res_coord, v_[1])
                        distances+=[dist]
                    dist_by_res[res] = distances   
                dist_by_pose["pose_{}".format(i)]= dist_by_res
            all_distances[k] = dist_by_pose
        # conf_pose_score keeps the total score of each pose
        min_distances = {}
        conf_pose_score = {}
        for k,v in all_distances.items():
            min_distances_by_pose = {}
            poses_score = []
            for k_, v_ in v.items():
                min_dist = {key: {'min_value': min(values), 'min_index': values.index(min(values))} for key, values in v_.items()}
                min_distances_by_pose[k_]=min_dist
                sum_min_dist = sum([min(v_) for v_ in v_.values()])
                poses_score += [sum_min_dist]
            min_distances[k] = min_distances_by_pose
            conf_pose_score[k]=poses_score
        return conf_pose_score
    
    def get_distance_score(self):
        conf_pose_score = self._get_distance_by_pose()
        all_scores = []
        conf_names = []
        for k,v in conf_pose_score.items():
            all_scores.append(v)
            for i in range(len(v)):
                name = "{}_pose{}".format(k, i)
                conf_names += [name]
        all_scores_flat = [i for i_ in all_scores for i in i_]
        mol_ids =[id.split("_")[0] for id in conf_names]
        df = pd.DataFrame({"id": mol_ids, "conf_id":conf_names, "distance_score": all_scores_flat})
        df["distance_rank"] = rankdata(df["distance_score"])
        return df

def main(output_folder):
    d = ProteinDocker(output_folder)
    df1 = d.get_docking_score()
    if os.path.exists(RESIDUE_COORDS):
        r = ResidueDistanceCalc(output_folder)
        df2 = r.get_distance_score()
        df = pd.merge(df1, df2, on=["id", "conf_id"], how="inner")
        df["rank"] = df.apply(lambda row: (row['distance_rank'] + row['docking_rank']) / 2, axis=1)
        df["idx_pose"] = df.index
        df["idx_pose"] = df.apply(lambda row: (row["idx_pose"]+1), axis=1)
        min_rank_indices = df.groupby('id')['rank'].idxmin()
        df_min = df.loc[min_rank_indices]
    else:
        df = df1
        df.rename(columns = {"docking_rank": "rank"}, inplace=True)
        df["idx_pose"] = df.index
        df["idx_pose"] = df.apply(lambda row: (row["idx_pose"]+1), axis=1)
        min_rank_indices = df.groupby('id')['rank'].idxmin()
        df_min = df.loc[min_rank_indices]
    df.to_csv(os.path.join(output_folder, "all_docking_results.csv"), index=False)
    df_min.to_csv(os.path.join(output_folder, "best_docking_results.csv"), index=False)

if __name__ == '__main__':
    output_folder = sys.argv[1]
    main(output_folder)

