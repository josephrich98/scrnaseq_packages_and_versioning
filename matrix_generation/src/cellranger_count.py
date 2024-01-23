import os
import sys
import yaml

import subprocess
import argparse
from datetime import datetime
import re
now = datetime.now()
date_directory_name = now.strftime('%y%m')

parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path) if parent_path not in sys.path else None

from scripts.organize_outputs import organize_output
from scripts.install_cellranger import install_cellranger_function

def ensure_correct_cellranger_version(instance):
    path_output_raw = subprocess.run("echo $PATH", shell=True, check=True, stdout=subprocess.PIPE, text=True)
    path_output = path_output_raw.stdout.strip()
    path_elements = path_output.split(':')
    new_path_elements = [element for element in path_elements if 'cellranger' not in element]

    desired_cellranger_path = f'{instance.package_path}/cellranger-{instance.cellranger_version}'
    
    # Define the regex pattern for the cellranger version
    current_cellranger_version_pattern = re.compile(r'cellranger-(\d+\.\d+\.\d+)')

    # Initialize variable to store the version
    current_cellranger_version = None
    all_cellranger_versions = []

    # Iterate through path elements to find the first cellranger version
    for element in path_elements:
        match = current_cellranger_version_pattern.search(element)
        if match:
            # Extract the version number
            if current_cellranger_version is None:
                current_cellranger_version = match.group(1)
            all_cellranger_versions.append(match.group(1))

    # Append the desired cellranger version at the beginning
    if current_cellranger_version != instance.cellranger_version:
        if instance.cellranger_version not in all_cellranger_versions:
            user_response = input(f"Cellranger {instance.cellranger_version} is not installed. Would you like to install it? (y/n): ")
            if user_response.lower() == 'y':
                install_cellranger_function(instance)
            else:
                print("Cellranger version in config.yaml not installed. To continue, please change config.yaml version or install this package")
                sys.exit("Installation aborted.")
        print(f"Successfully added cellranger version {instance.cellranger_version} as the current version in PATH.")
    else:
        print(f"Cellranger version {instance.cellranger_version} is correct")

    # Join the elements back into a PATH string
    new_path_elements.insert(0, desired_cellranger_path)
    updated_path = ':'.join(new_path_elements)
    os.environ["PATH"] = updated_path

def cellranger_install_transcriptome_function(instance):
    cellranger_reference_path = os.path.join(instance.reference_path, "cellranger", date_directory_name)
    if not os.path.exists(cellranger_reference_path):
        os.makedirs(cellranger_reference_path)
    os.chdir(cellranger_reference_path)
    if not os.listdir(cellranger_reference_path):
        subprocess.run(f"wget {instance.cellranger_transcriptome_link}", shell=True, executable="/bin/bash")
        subprocess.run(f"tar -xzvf {instance.cellranger_transcriptome_link.split('/')[-1]}")
        

def find_transcriptome(instance, reference_path):
    numbers = []

    # 1. Find the most recent transcriptome
    for file in os.listdir(reference_path):
        try:
            num = int(file)
            numbers.append(num)
        except ValueError:
            continue

    # If we found transcriptomes
    if numbers:
        if instance.reference_selector == "most_recent":
            transcriptome_date = max(numbers)
        else:
            transcriptome_date = instance.reference_selector
            
        grch38_dir = os.path.join(
            reference_path, str(transcriptome_date), "grch38_transcriptome"
        )
        
        # 2. Look within the grch38_transcriptome of the most recent transcriptome
        if os.path.exists(grch38_dir):
            for child_dir in os.listdir(grch38_dir):
                # 3. Check for directories starting with refdata
                if child_dir.startswith("refdata"):
                    return os.path.join(grch38_dir, child_dir)

    # Return None if we didn't find any matching directory
    return None


def cellranger_count_function(instance, baseline):
    if baseline:
        seed_list = (0,)
        frac_list = (1.0,)
    else:
        seed_list = instance.seed_list
        frac_list = instance.frac_list

    if instance.cellranger_version != "":
        ensure_correct_cellranger_version(instance)

    cellranger_version_str = str(instance.cellranger_version).replace('.', '_')

    cellranger_output_directory = os.path.join(instance.output_directory, f"cellranger{cellranger_version_str}")
    os.makedirs(cellranger_output_directory, exist_ok=True)
    os.chdir(cellranger_output_directory)
    cellranger_reference_path = os.path.join(
        os.path.dirname(instance.main_directory),
        "reference_files",
        "cellranger",
    )
    
    transcriptome_path = find_transcriptome(instance, cellranger_reference_path)

    for frac in frac_list:
        for seed in seed_list:
            frac_str = str(frac).replace('.', '_')
            
            if baseline or int(frac) == 1:
                specific_fastq_directory = instance.original_fastq_directory
            else:
                specific_fastq_directory = os.path.join(instance.downsampled_fastq_directory, f"frac{frac_str}_seed{seed}")
            
            command = [
                "cellranger",
                "count",
                f"--id=frac{frac_str}_seed{seed}",
                f"--transcriptome={transcriptome_path}",
                f"--fastqs={specific_fastq_directory}",
                f"--sample={instance.data_name}",
                "--localcores=8",
                "--localmem=64",
            ]  # modify ID to change output folder name

            subprocess.run(command)

            print(f"cellranger count done for seed {seed} and frac {frac}")
            
            organize_output(instance.output_directory, seed, frac_str, matrix_source="cellranger", matrix_version = cellranger_version_str)

    print(f"cellranger count done!")


if __name__ == "__main__":    
    from fastq_processor import FastqProcessor
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--baseline', action='store_true', help='True if using baseline data (will override any seeds and counts accordingly), False if using downsampled fastq data')
    args = parser.parse_args()
   
    with open(f'{parent_path}/config.yaml', 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)

    fastq_processor.cellranger_count(baseline=args.baseline)
