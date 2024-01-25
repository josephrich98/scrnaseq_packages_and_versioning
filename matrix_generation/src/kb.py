import os
import sys
import yaml
import pkg_resources

parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_path) if parent_path not in sys.path else None

import subprocess
import json
import shutil
from datetime import datetime
import argparse
from scripts.organize_outputs import organize_output
now = datetime.now()


def check_kb_version(target_version):
    version = pkg_resources.get_distribution("kb-python").version
    if version == target_version:
        print(f"kb version {version} matches config.yaml. Proceeding to run kb count")
    else:
        raise ValueError(f"kb target version = {target_version}, but actual kb version = {version}. Please check config.yaml and/or the python environment. kb-python can be installed with pip install kb-python gget ffq")

def custom_sort(filename):
    # Define order for file types
    file_type_order = {'R1': 0, 'R2': 1, 'I1': 2, 'I2': 3}

    # Split filename by '_' to extract file type and lane information
    parts = filename.split("_")

    # Extract lane number; assuming lane info is of the format 'L00X'
    lane = int(parts[-3][1:4])  # e.g., extracts '001' from 'L001'

    # Get the order value for the file type, e.g., 'R1'
    file_type = parts[-2].split(".")[0]  # e.g., extracts 'R1' from 'R1_001.fastq.gz'

    # Return a tuple: (file_type_order, lane_number)
    return (lane, file_type_order.get(file_type, 999))  # 999 as a default value for any unexpected file types




def find_kb_reference(instance, count = False):
    kb_reference_path_parent = os.path.join(instance.reference_path, "kb")
    if instance.reference_selector == "most_recent":
        if count:
            directories = [item for item in os.listdir(kb_reference_path_parent) if os.path.isdir(os.path.join(kb_reference_path_parent, item))]
            kb_reference_path = max(int(dir_name.split('_')[0]) for dir_name in directories if '_' in dir_name)

        else:
            kb_reference_path = os.path.join(kb_reference_path_parent, now.strftime('%y%m'))
    else:
        kb_reference_path = os.path.join(kb_reference_path_parent, instance.reference_selector)
    
    kb_version_parts = instance.kb_version.split('.')
    
    kb_reference_path_full = os.path.join(kb_reference_path_parent, str(kb_reference_path) + f"_v{kb_version_parts[1]}")
    
    if count:
        if not os.path.exists(kb_reference_path_full):
            raise Exception(f"{kb_reference_path_full} does not exist. Please run kb ref with this date and version, or check values in config.yaml")

    return kb_reference_path_full

def kb_ref_function(instance):
    ensembl_path = os.path.join(instance.reference_path, "ensembl", str(instance.ensembl_release))
    
    if not os.path.exists(ensembl_path):
        os.makedirs(ensembl_path)
    os.chdir(ensembl_path)
    if not os.listdir(ensembl_path):
        gget_command = f"gget ref -o {ensembl_path}/info.json -r {instance.ensembl_release} -d -w dna,gtf homo_sapiens"
        print(gget_command)
        subprocess.run(gget_command, shell=True, executable="/bin/bash")
    
    kb_reference_path = find_kb_reference(instance)
    if not os.path.exists(kb_reference_path):
        os.makedirs(kb_reference_path)
    os.chdir(kb_reference_path)
    if instance.kb_version != "":
        check_kb_version(instance.kb_version)
    
    if instance.kb_workflow == "nac":
        kb_ref_command = f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/fasta_spliced.fasta -f2 {kb_reference_path}/fasta_unspliced.fasta -c1 {kb_reference_path}/c_spliced.txt -c2 {kb_reference_path}/c_unspliced.txt --workflow=nac {ensembl_path}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz {ensembl_path}/Homo_sapiens.GRCh38.{instance.ensembl_release}.gtf.gz"
    else:   
        kb_ref_command = f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/transcriptome.fa {ensembl_path}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz {ensembl_path}/Homo_sapiens.GRCh38.{instance.ensembl_release}.gtf.gz"
    
    print(kb_ref_command)
    subprocess.run(kb_ref_command, shell=True, executable="/bin/bash")
    
        
def kb_count_function(instance, baseline, threads):
    if baseline:
        seed_list = (0,)
        frac_list = (1.0,)
    else:
        seed_list = instance.seed_list
        frac_list = instance.frac_list

    if instance.kb_version != "":
        check_kb_version(instance.kb_version)

    print("about to enter kb count")
    for frac in frac_list:
        for seed in seed_list:
            frac_str = str(frac).replace('.', '_')
            
            if baseline or int(frac) == 1:
                specific_fastq_directory = instance.original_fastq_directory
            else:
                specific_fastq_directory = os.path.join(instance.downsampled_fastq_directory, f"frac{frac_str}_seed{seed}")

            kb_version_str = str(instance.kb_version).replace('.', '_')

            out_dir = os.path.join(instance.output_directory, f"kb{kb_version_str}", f"frac{frac_str}_seed{seed}")
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            # else:
            #     response = input(f"Warning: Output folder for frac{frac_str}_seed{seed} already exists in this project. Do you want to continue? (y/n) ")
            #     if response.lower() != "y":
            #         print("Operation aborted by the user.")
            #         sys.exit()
            
            os.chdir(specific_fastq_directory)

            print(f"Current directory: {specific_fastq_directory}")
            
            kb_reference_path = find_kb_reference(instance, count = True)

            # Build the command
            kb_count_command = [
                'kb', 'count', 
                # '--filter',
                # '--gene-names',
                '-i', f'{kb_reference_path}/index.idx', 
                '-g', f'{kb_reference_path}/t2g.txt', 
                '-x', f'{instance.kb_sequencing_technology}', 
                '-o', f"{out_dir}",
                '-t', str(threads)
            ]

            if instance.kb_count_format == "h5ad":
                kb_count_command.append("--h5ad")

            if instance.kb_workflow == "nac":
                # kb_count_command.append(f"--workflow=nac -c1 {kb_reference_path}/c_spliced.txt -c2 {kb_reference_path}/c_unspliced.txt --sum=total")
                new_arguments = ["--workflow=nac", "-c1", f"{kb_reference_path}/c_spliced.txt", "-c2", f"{kb_reference_path}/c_unspliced.txt", "--sum=total"]
                for argument in new_arguments:
                    kb_count_command.append(argument)
            
            # Add the fastq files
            fastq_extensions = ('.fq', '.fastq', '.fq.gz', '.fastq.gz')
            fastq_files = [filename for filename in os.listdir(specific_fastq_directory) if filename.endswith(fastq_extensions)]
            
            sorted_files = sorted(fastq_files, key=custom_sort)

            sorted_files = [file for file in sorted_files if "_I1_" not in file and "_I2_" not in file]

            kb_version = pkg_resources.get_distribution("kb-python").version

            kb_version_major = int(kb_version.split('.')[1])

            if kb_version_major >= 28:
                lane_files = {}
                for file in sorted_files:
                    lane_number = file.split('_')[-3]  # Extract lane number (e.g., L002)
                    lane_files.setdefault(lane_number, []).append(file)

                with open('batch.txt', 'w') as f:
                    for lane, files in lane_files.items():
                        line = lane + "\t" + '\t'.join(files) + "\n"
                        f.write(line)

                kb_count_command.append("--batch-barcodes batch.txt")

            else:
                # Append sorted filenames to kb_count_command
                for file in sorted_files:
                    kb_count_command.append(file)
            
            print(' '.join(kb_count_command))

            kb_count_command = ' '.join(kb_count_command)

            # Run the command
            subprocess.run(kb_count_command, shell=True, executable="/bin/bash")

            # organize_output(instance.output_directory, seed, frac_str, matrix_source = "kb")
            
            print(f"kb count finished for seed {seed} and frac {frac_str}")
    
    for frac in frac_list:
        for seed in seed_list:
            frac_str = str(frac).replace('.', '_')
            organize_output(instance.output_directory, seed, frac_str, matrix_source = "kb", matrix_version = kb_version_str)

    print("kb count finished")




    
            
if __name__ == "__main__":
    from fastq_processor import FastqProcessor
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--baseline', action='store_true', help='True if using baseline data (will override any seeds and counts accordingly), False if using downsampled fastq data')
    parser.add_argument('-t', '--threads', default='8', help='Number of threads for kb count')
    parser.add_argument('-r', '--run_ref', action='store_true', help='Run kb ref in addition to kb count')
    args = parser.parse_args()
   
    with open(f'{parent_path}/config.yaml', 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)

    if args.run_ref:
        fastq_processor.kb_ref()
    
    fastq_processor.kb_count(baseline=args.baseline, threads=args.threads)
