import os
import sys
import shutil
import subprocess
import yaml

parent_path = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))
)  # gets the parent directory of the file that this script is in
sys.path.append(parent_path)  # adds the parent directory to the path

def install_cellranger_function(instance):
    subprocess.run(f"""
        cd {instance.package_path} &&
        
        echo 'curl -o cellranger-{instance.cellranger_version}.tar.gz "{instance.cellranger_package_link}"' &&
        curl -o cellranger-{instance.cellranger_version}.tar.gz "{instance.cellranger_package_link}" &&
        
        echo 'tar -xzvf cellranger-{instance.cellranger_version}.tar.gz' &&
        tar -xzvf cellranger-{instance.cellranger_version}.tar.gz &&
        
        export PATH={instance.package_path}/cellranger-{instance.cellranger_version}:$PATH
    """, shell=True, executable="/bin/bash")
    with open(f"{os.path.expanduser('~')}/.bashrc", "a") as f:
        f.write(f"\nexport PATH={instance.package_path}/cellranger-{instance.cellranger_version}:$PATH")

if __name__ == "__main__":
    from fastq_processor import FastqProcessor
   
    with open(f'{parent_path}/config.yaml', 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)
    
    fastq_processor.install_cellranger()