# SenSaaS-flex

This package is an adaptation from the original [SenSaaS](https://github.com/SENSAAS/sensaas) package, described in [Biyuzan et al, 2024](https://academic.oup.com/bioinformatics/article/40/3/btae105/7612231). We have downloaded the [code](https://chemoinfo.ipmc.cnrs.fr/Sensaas-flex/sensaas-flex-main.tar.gz) on 1st July 2024 and adapted it for our use. If you use it, please cite the original authors.

## Installation

This installation instructions are only tested in Ubuntu 22.04 LTS. 

'''bash
conda create -n sensaas python=3.8
conda activate sensaas
git clone https://github.com/ersilia-os/sensaas-flex
cd sensaas-flex
wget https://anaconda.org/open3d-admin/open3d/0.12.0/download/linux-64/open3d-0.12.0-py38_0.tar.bz2
conda install open3d-0.12.0-py38_0.tar.bz2
pip install -r requirements.txt
'''

## Usage
To align the conformers of a database of molecules to one query molecule (target) that stays fixed, simply run:
'''bash
python meta-sensaasflex.py [query.sdf] [database.sdf] [output_folder]
'''
The molecules will need to be passed as conformers in an sdf file (see: https://github.com/ersilia-os/smiles-to-3d)

## License
This package is licensed under a BSD 3-Clause