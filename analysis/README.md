# This is the repository to compare Seurat and Scanpy, count matrix generation methods, and full sized and downsampled matrices

# To run docker container:
docker run --name seurat_vs_scanpy -d -p 8787:8787 -e PASSWORD=yourpassword -v /path/to/code:/workspace josephrich98/seurat_vs_scanpy_analysis
- visit http://localhost:8787
- sign in with username rstudio, password yourpassword (set above)

Activate the appropriate conda environment
analysis_env: Scanpy 1.9.5
sc14_umap5: Scanpy 1.4.6, UMAP-learn 0.5.1
sc14_umap4: Scanpy 1.4.6, UMAP-learn 0.4.6


# To run in conda environment: build the appropriate conda environment from yaml file
Installations needed to be performed manually (cannot be stored in yaml due to dependency conflicts or difficulty finding package upon building):
all environments: pip install git+https://github.com/has2k1/scikit-misc.git@269f61e
sc14_umap4: pip install kb-python==0.27.3 ngs-tools==1.8.5 --no-deps