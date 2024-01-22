import os
import shutil

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