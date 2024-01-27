import os
import sys
import shutil
import argparse

def organize_output(src_dir, seed, frac_str, matrix_source, matrix_version = ""):
    if matrix_source == "cellranger":
        full_src_dir = os.path.join(src_dir, f"{matrix_source}{matrix_version}", f"frac{frac_str}_seed{seed}", "outs", "raw_feature_bc_matrix")
        # full_src_dir_filtered = os.path.join(src_dir, matrix_source, f"frac{frac_str}_seed{seed}", "outs", "filtered_feature_bc_matrix")
    elif matrix_source == "kb":
        full_src_dir = os.path.join(src_dir, f"{matrix_source}{matrix_version}", f"frac{frac_str}_seed{seed}", "counts_unfiltered")
    dest_dir_root = os.path.join(src_dir, f"{matrix_source}{matrix_version}_raw_feature_bc_matrix_collection")
    if not os.path.exists(dest_dir_root):
        os.makedirs(dest_dir_root, exist_ok=True)
        
    dest_dir = os.path.join(dest_dir_root, f"frac{frac_str}_seed{seed}")
    
    shutil.copytree(full_src_dir, dest_dir)

if __name__ == "__main__":
    parent_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(parent_path) if parent_path not in sys.path else None
    
    from fastq_processor import FastqProcessor
    import yaml

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--baseline', action='store_true', help='Whether or not to include full data in fraction list')
    parser.add_argument('-v', '--version', default = "", help='Provide version of kb/cellranger used')
    parser.add_argument('-t', '--tool', choices=['kb', 'cellranger'], default = 'kb', help='The tool to use, either "kb" or "cellranger"')
    args = parser.parse_args()
   
    with open(f'{parent_path}/config.yaml', 'r') as file:
        config = yaml.safe_load(file)
        
    fastq_processor = FastqProcessor(**config)

    if args.baseline:
        fastq_processor.frac_list.append(1.0)

    if args.version:
        matrix_version = args.version.replace('.', '_')

    for seed in fastq_processor.seed_list:
        for frac in fastq_processor.frac_list:
            frac_str = str(frac).replace('.', '_')
            if frac_str == "1_0":
                seed = 0
            organize_output(fastq_processor.output_directory, seed, frac_str, matrix_source = args.tool, matrix_version = matrix_version)
