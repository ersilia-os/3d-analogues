import os
import argparse
import subprocess

root = os.path.dirname(os.path.abspath(__file__))

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
            generate_conformers(query_input, os.path.join(output_folder, "conformer_query.sdf"))
        else:
            from conformers import ObabelConf
            obabelconf = ObabelConf()
            obabelconf.generate_conformers(sampled_input, os.path.join(output_folder, "conformers.sdf"))
            obabelconf.generate_conformers(query_input, os.path.join(output_folder, "conformer_query.sdf"))
    
        cmd = f"python {root}/docking.py {output_folder}"
        subprocess.Popen(cmd, shell=True).wait()


    

def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Evaluates analogues in the 3D space')

    parser.add_argument('-q',
                        dest='query_file',
                        required=False,
                        metavar='<file>',
                        help='Query molecule in CSV format. Must have a smiles and id column')
    
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