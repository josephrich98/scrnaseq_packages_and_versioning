# This is the repository to compare Seurat and Scanpy, count matrix generation methods, and full sized and downsampled matrices

# To run docker container:
docker run --name seurat_vs_scanpy -d -p 8787:8787 -e PASSWORD=yourpassword -v /path/to/code:/workspace josephrich98/seurat_vs_scanpy_analysis
- visit http://localhost:8787
- sign in with username rstudio, password yourpassword (set above)

conda activate seurat_vs_scanpy_analysis