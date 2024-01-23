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
    if instance.reference_selector == "most_recent":
        if count:
            directories = [d for d in os.listdir(instance.reference_path) if os.path.isdir(os.path.join(instance.reference_path, d))]
            kb_reference_path = max(int(dir_name.split('_')[0]) for dir_name in directories if '_' in dir_name)

        else:
            kb_reference_path = os.path.join(instance.reference_path, "kb", now.strftime('%y%m'))
    else:
        kb_reference_path = os.path.join(instance.reference_path, "kb", instance.reference_selector)
    
    kb_version_parts = instance.kb_version.split('.')

    kb_reference_path_full = kb_reference_path + f"_v{kb_version_parts[1]}"

    if count:
        if not os.path.exists(kb_reference_path_full):
            raise Exception(f"{kb_reference_path_full} does not exist. Please run kb ref with this date and version, or check values in config.yaml")
    
    return kb_reference_path_full

def kb_ref_function(instance):
    kb_reference_path = find_kb_reference(instance)
    if not os.path.exists(kb_reference_path):
        os.makedirs(kb_reference_path)
    os.chdir(kb_reference_path)
    if instance.kb_version != "":
        check_kb_version(instance.kb_version)
    
    kb_major_version = instance.kb_version.split('.')[1]
    
    if int(kb_major_version) < 27:
        if instance.reference_selector == "most_recent":
            subprocess.run(f"wget $(gget ref --ftp -w dna homo_sapiens)", shell=True, executable="/bin/bash") 
            subprocess.run(f"wget $(gget ref --ftp -w gtf homo_sapiens)", shell=True, executable="/bin/bash") 
        else:
            subprocess.run(f"wget $(gget ref -r {instance.ensembl_release} --ftp -w dna homo_sapiens)", shell=True, executable="/bin/bash")
            subprocess.run(f"wget $(gget ref -r {instance.ensembl_release} --ftp -w gtf homo_sapiens)", shell=True, executable="/bin/bash") 
        kb_ref_command = f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/transcriptome.fa {kb_reference_path}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz {kb_reference_path}/Homo_sapiens.GRCh38.110.gtf.gz"
    else:
        if instance.kb_workflow == "nac":
            if instance.reference_selector == "most_recent":
                kb_ref_command = f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/fasta_spliced.fasta -f2 {kb_reference_path}/fasta_unspliced.fasta -c1 {kb_reference_path}/c_spliced.txt -c2 {kb_reference_path}/c_unspliced.txt --workflow=nac $(gget ref --ftp -w dna,gtf homo_sapiens)"
            else:
                # subprocess.run(f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/fasta_spliced.fasta -f2 {kb_reference_path}/fasta_unspliced.fasta -c1 {kb_reference_path}/c_spliced.txt -c2 {kb_reference_path}/c_unspliced.txt --workflow=nac {instance.kb_fasta_link} {instance.kb_gtf_link}", shell=True, executable="/bin/bash")
                kb_ref_command = f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/fasta_spliced.fasta -f2 {kb_reference_path}/fasta_unspliced.fasta -c1 {kb_reference_path}/c_spliced.txt -c2 {kb_reference_path}/c_unspliced.txt --workflow=nac $(gget ref --ftp -r {instance.ensembl_release} -w dna,gtf homo_sapiens)"
        else:
            if instance.reference_selector == "most_recent":
                kb_ref_command = f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/transcriptome.fa $(gget ref --ftp -w dna,gtf homo_sapiens)"
            else:
                # subprocess.run(f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/transcriptome.fa {instance.kb_fasta_link} {instance.kb_gtf_link}", shell=True, executable="/bin/bash")
                kb_ref_command = f"kb ref -i {kb_reference_path}/index.idx -g {kb_reference_path}/t2g.txt -f1 {kb_reference_path}/transcriptome.fa $(gget ref -r {instance.ensembl_release} --ftp -w dna,gtf homo_sapiens)"
    
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
                '-t', threads
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

            # Append sorted filenames to kb_count_command
            for file in sorted_files:
                kb_count_command.append(file)
            
                    
            print(' '.join(kb_count_command))
            # Run the command
            result = subprocess.run(kb_count_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Check for errors
            if result.returncode != 0:
                print(f"Command failed with error: {result.stderr.decode()}")
            else:
                print(result.stdout.decode())
                
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
