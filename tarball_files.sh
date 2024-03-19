#!/bin/bash

# Change to the specified directory
cd "/Users/joeyrich/Desktop/local/seurat_scanpy_final/scrnaseq_packages_and_versioning/analysis/count_matrix_collection/SC3_v3_NextGem_SI_PBMC_10K" || exit

# Loop through each directory in the current directory
for parent_dir in */ ; do
    # Trim the trailing slash
    parent_dir=${parent_dir%/}
    # Enter the directory
    cd "$parent_dir"

    # Loop through each subdirectory
    for sub_dir in */ ; do
        # Trim the trailing slash
        sub_dir=${sub_dir%/}
        # Hold the directory name
        new_name="${parent_dir}_${sub_dir}"
        # Rename the folder
        mv "$sub_dir" "$new_name"
        # Tarball the folder
        tar -czvf "${new_name}.tar.gz" --exclude='.DS_Store' "$new_name"
        # Rename the subdirectory back
        mv "$new_name" "$sub_dir"
    done

    # Go back to the parent directory
    cd ..
done
