import os
import sys
import argparse
import subprocess
import shutil
import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
abspath = os.path.abspath(__file__)
sys.path.append(abspath)

import docking
from utils import extract_best_conformer

def main() -> None:
    args = parseArgs()

    query_input = args.query_file
    sampled_input = args.sample_file
    output_folder = args.output_folder
    cdpkit = args.cdpkit

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
        print("Output Folder not existing, it will be created")

    if query_input is None:
        if cdpkit=="True":
            print("HERE")
            from smiles3d.smiles3d import generate_conformers
            generate_conformers(sampled_input, os.path.join(output_folder, "conformers.sdf"))
        else:
            from conformers import ObabelConf
            obabelconf = ObabelConf()
            obabelconf.generate_conformers(sampled_input, os.path.join(output_folder, "conformers.sdf"))
    
        cmd = f"python {root}/docking.py {output_folder}"
        subprocess.Popen(cmd, shell=True).wait()

    else: 
        if cdpkit=="True":
            from smiles3d.smiles3d import generate_conformers
            generate_conformers(sampled_input, os.path.join(output_folder, "conformers.sdf"))

            if query_input.endswith(".sdf"):
                output_file = os.path.join(output_folder, "conformer_query.sdf")
                shutil.copy(query_input, output_file)
            else:
                generate_conformers(query_input, os.path.join(output_folder, "conformer_query.sdf"))
        else:
            from conformers import ObabelConf
            obabelconf = ObabelConf()
            obabelconf.generate_conformers(sampled_input, os.path.join(output_folder, "conformers.sdf"))

            if query_input.endswith(".sdf"):
                output_file = os.path.join(output_folder, "conformer_query.sdf")
                shutil.copy(query_input, output_file)
            else:
                obabelconf.generate_conformers(query_input, os.path.join(output_folder, "conformer_query.sdf"))
        

        #we want to keep only the active conformer as a query - if not inputed directly, calculate with docking
        docking.run(output_folder, "conformer_query.sdf")
        #Now we do docking with all the sampled smiles
        docking.run(output_folder, "conformers.sdf")

        best_conf_query = pd.read_csv(os.path.join(output_folder, "best_docking_results_query.csv"))
        best_conf_query = best_conf_query["conf_id"].tolist()
        extract_best_conformer(os.path.join(output_folder, "conformer_query.sdf"), best_conf_query, os.path.join(output_folder, "conformer_query_best.sdf"))
        best_conf = pd.read_csv(os.path.join(output_folder, "best_docking_results.csv"))
        best_conf = best_conf["conf_id"].tolist()
        extract_best_conformer(os.path.join(output_folder,"conformers.sdf"), best_conf, os.path.join(output_folder,"conformers_best.sdf"))


        from vsflow3d.vsflow3d import run_vsflow
        run_vsflow(os.path.join(output_folder, "conformer_query_best.sdf"), os.path.join(output_folder, "conformers_best.sdf"), output_folder)

        # rename vsflow results for easier merging
        vsflow_results = pd.read_csv(os.path.join(output_folder, "vsflow_results.csv"))
        ids = [id.split("_")[0] for id in vsflow_results["id"]]
        vsflow_results.rename(columns={"id":"conf_id"}, inplace=True)
        vsflow_results["id"] = ids
        docking_results = pd.read_csv(os.path.join(output_folder, "best_docking_results.csv"))
        df = pd.merge(docking_results, vsflow_results, on=["id", "conf_id"], how = "left")
        df = df[['id', 'conf_id', 'pose_id', 'smiles', 'querysmiles','docking_score', 'docking_rank', 'distance_score',
                 'distance_rank', 'rank', 'idx_pose', 'ComboScore', 'ShapeSim', 'Fp3dSim']]
        df.to_csv(os.path.join(output_folder, "all_results.csv"), index=False)
    

def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Evaluates analogues in the 3D space')

    parser.add_argument('-q',
                        dest='query_file',
                        required=False,
                        metavar='<file>',
                        help='Query molecule in CSV format. Must have a smiles and id column. Can also be a .sdf file with a single conformer.')
    
    parser.add_argument('-s',
                        dest='sample_file',
                        required=True,
                        metavar='<file>',
                        help='Sampled molecules in CSV format. Must have a smiles and id column')
    
    parser.add_argument('-o',
                        dest='output_folder',
                        required=True,
                        metavar='<folder>',
                        help='Output folder.')
    
    parser.add_argument('-cdpkit',
                        dest='cdpkit',
                        required=False,
                        default=False,
                        help='Use the CPDKit to generate 3D conformers.')
    
    return parser.parse_args()


if __name__ == '__main__':
    main()