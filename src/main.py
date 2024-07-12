import os
import sys
import argparse
import subprocess
import shutil
import pandas as pd
from smiles3d.smiles3d import generate_conformers

root = os.path.dirname(os.path.abspath(__file__))
abspath = os.path.abspath(__file__)
sys.path.append(abspath)

import docking
from utils import SDFProcessor
from vsflow3d import run_vsflow

def main() -> None:
    args = parseArgs()

    query_input = args.query_file
    sampled_input = args.sample_file
    output_folder = args.output_folder
    cdpkit = args.cdpkit

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
        print("Output Folder not existing, it will be created")
        
    #STEP 1: Generate Conformers
    if cdpkit=="True":
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
    
    #STEP 2: Select best query conformer if not passed.
    #if an SDF file with only one conformer is available, this will run for that molecule 
    docking.run(output_folder, "conformer_query.sdf")
    best_conf_query = pd.read_csv(os.path.join(output_folder, "best_docking_results_query.csv"))
    best_conf_query = best_conf_query["pose_id"].tolist()
    sp = SDFProcessor(output_folder)
    sp.add_hs("docking_result_query.sdf","docking_result_query_hs.sdf")
    sp.add_pose_to_name("docking_result_query_hs.sdf", "docking_result_pose_query.sdf")
    best_conf_query = pd.read_csv(os.path.join(output_folder, "best_docking_results_query.csv"))
    best_conf_query = best_conf_query["pose_id"].tolist()
    sp.extract_best_conformer("docking_result_pose_query.sdf", "docking_query_best.sdf", best_conf_query)
    os.remove(os.path.join(output_folder, "docking_result_query.sdf"))
    os.remove(os.path.join(output_folder, "docking_result_query_hs.sdf"))

    #STEP 3: Align Conformers from Sample Molecules to Query
    cmd = f"python {os.path.join(root, 'sensaasflex', 'meta-sensaasflex.py')} {os.path.join(output_folder,'docking_query_best.sdf')} {os.path.join(output_folder,'conformers.sdf')} {output_folder}"
    subprocess.Popen(cmd, shell=True).wait()
    original_names = sp.get_original_names("conformers.sdf")
    df = sp.create_matrix_csv(original_names)
    best_ids = df["conf_id"].tolist()
    sp.append_names_to_sdf("catsensaas.sdf", "catsensaas_names.sdf", original_names)
    sp.extract_best_conformer("catsensaas_names.sdf", "catsensaasbest.sdf", best_ids)
    os.remove(os.path.join(output_folder, "matrix-sensaas.txt"))
    os.remove(os.path.join(output_folder, "catsensaas.sdf"))
    os.remove(os.path.join(output_folder, "bestsensaas.sdf"))
    #optional: add a filter for Sensaas Score minimum

    #STEP 4: Docking
    docking.run(output_folder, "catsensaasbest.sdf")
    df = pd.read_csv(os.path.join(output_folder, "all_docking_results.csv"))
    ids = df["pose_id"].tolist()
    sp.append_names_to_sdf('docking_result.sdf', 'docking_result_names.sdf',ids)

    best_df = pd.read_csv(os.path.join(output_folder, "best_docking_results.csv"))
    best_ids = best_df["pose_id"].tolist()
    sp.extract_best_conformer('docking_result_names.sdf', "docking_result_best.sdf", best_ids)
    os.remove(os.path.join(output_folder, "docking_result.sdf"))

    #STEP 5: perform 3D Similarity using VSFlow
    run_vsflow(os.path.join(output_folder, "docking_query_best.sdf"), os.path.join(output_folder, "catsensaasbest.sdf"), output_folder)

    #STEP6: Join all results into single file
    # rename vsflow results for easier merging
    vsflow_results = pd.read_csv(os.path.join(output_folder, "vsflow_results.csv"))
    ids = [id.split("_")[0] for id in vsflow_results["id"]]
    vsflow_results.rename(columns={"id":"conf_id"}, inplace=True)
    vsflow_results["id"] = ids
    docking_results = pd.read_csv(os.path.join(output_folder, "best_docking_results.csv"))
    df = pd.merge(docking_results, vsflow_results, on=["id", "conf_id"], how = "left")
    sensaas_results = pd.read_csv(os.path.join(output_folder, "matrix_sensaas.csv"))
    sensaas_results= sensaas_results.drop(columns=["idx_pose"])
    df_ = pd.merge(df, sensaas_results, on = ["id", "conf_id"], how="left")
    df_ = df_[['id', 'conf_id', 'pose_id', 'smiles', 'querysmiles','sensaas_score','docking_score', 'docking_rank', 'distance_score',
                'distance_rank', 'rank', 'ComboScore', 'ShapeSim', 'Fp3dSim']]
    df_.to_csv(os.path.join(output_folder, "all_results.csv"), index=False)


def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Evaluates analogues in the 3D space')

    parser.add_argument('-q',
                        dest='query_file',
                        required=False,
                        metavar='<file>',
                        help='Query molecule in CSV format. Must have a smiles and id column. Can also be a .sdf file with a single conformer.')
    
    parser.add_argument('-s',
                        dest='sample_file',
                        required=False,
                        metavar='<file>',
                        help='Sampled molecules in CSV format. Must have a smiles and id column')
    
    parser.add_argument('-o',
                        dest='output_folder',
                        required=True,
                        metavar='<folder>',
                        help='Output folder.')
    
    parser.add_argument('--cdpkit',
                        dest='cdpkit',
                        required=False,
                        default=True,
                        help='Use the CPDKit to generate 3D conformers.')
    
    return parser.parse_args()

if __name__ == '__main__':
    main()