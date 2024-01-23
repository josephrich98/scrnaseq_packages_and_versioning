import os
import sys
import yaml

parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path) if parent_path not in sys.path else None

from scripts.install_cellranger import install_cellranger_function

import shutil
import argparse
import subprocess
from datetime import datetime
now = datetime.now()
date_directory_name = now.strftime('%y%m')
    
def project_setup_function(instance):    
    if not os.path.exists(instance.package_path):
        os.makedirs(instance.package_path)

    os.chdir(instance.package_path)

    ### PACKAGE DOWNLOADS ###
    # Seqtk
    if shutil.which("seqtk") is None:
        subprocess.run(f"""
            cd {instance.package_path} &&
            git clone https://github.com/lh3/seqtk.git &&
            cd seqtk &&
            make &&
            chmod +x {instance.package_path}/seqtk/seqtk
            export PATH=$PATH:{instance.package_path}/seqtk
        """, shell=True, executable="/bin/bash")
        os.environ["PATH"] = f":{instance.package_path}/seqtk" + os.environ["PATH"]
        with open(f"{os.path.expanduser('~')}/.bashrc", "a") as f:
            f.write(f"\nexport PATH={instance.package_path}/seqtk:$PATH")

    # cellranger
    if shutil.which("cellranger") is None and instance.count_matrix_generation_method != "kb":
        install_cellranger_function(instance)

    ### Directory organization ###
    if not os.path.exists(instance.reference_path):
        os.makedirs(instance.reference_path)
        os.makedirs(f"{instance.reference_path}/kb/{date_directory_name}")
        os.makedirs(f"{instance.reference_path}/cellranger/{date_directory_name}")

    if not os.path.exists(instance.output_directory):
        os.makedirs(instance.output_directory)
        os.makedirs(f"{instance.output_directory}/kb")
        os.makedirs(f"{instance.output_directory}/cellranger")
        os.makedirs(f"{instance.output_directory}/10x_genomics")

    if not os.path.exists(instance.data_directory):
        os.makedirs(instance.downsampled_fastq_directory)
        os.makedirs(instance.original_fastq_directory)

    print("Setup finished. Please source your bashrc file or start a new terminal for seqtk and cellranger commands to take effect")
    
if __name__ == "__main__":
    parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(parent_path) if parent_path not in sys.path else None
    
    from fastq_processor import FastqProcessor
   
    with open(f'{parent_path}/config.yaml', 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)
    
    fastq_processor.project_setup()

