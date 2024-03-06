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

To perform analysis, select a yaml file (with modifications as desired) and run the appropriate Rmd notebook. Each notebook is also present as a Jupyter notebook (not compatible with the docker image), which can be run locally or in Google Colab (links below). A yaml file is present for each figure. The appropriate notebook to which each yaml files belongs is as follows:

- aggregate_plots: Fig2_Supp_Fig11_Supp_Fig12
- Seurat_v_Scanpy: Fig1, Supp_Fig2, Supp_Fig3, Supp_Fig4, Supp_Fig5
- Seurat_v_Seurat: Supp_Fig7, Supp_Fig9, Supp_Fig13, Supp_Fig14, Fig3_cellranger, Fig3_kallisto
- Scanpy_v_Scanpy: Supp_Fig8, Supp_Fig10
- Scanpy_version_comparison: Fig3_scanpy
- Seurat_version_comparison: Fig3_seurat

Note: Supp Fig 6 (Extended UMAP analysis) is derived from portions of the runs of Fig1, Supp_Fig2, Supp_Fig3, and Supp_Fig4. Supp Fig 15 (UMI filtering analysis) is derived from portions of the runs of Supp_Fig7, Supp_Fig13, and Supp_Fig14.


Google colab links:

aggregate_plots: [![Run in Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/josephrich98/scrnaseq_packages_and_versioning/blob/main/analysis/ipynb/aggregate_plots.ipynb)

Scanpy_v_Scanpy: [![Run in Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/josephrich98/scrnaseq_packages_and_versioning/blob/main/analysis/ipynb/Scanpy_v_Scanpy.ipynb)

Scanpy_version_comparison: [![Run in Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/josephrich98/scrnaseq_packages_and_versioning/blob/main/analysis/ipynb/Scanpy_version_comparison.ipynb)

Seurat_v_Scanpy: [![Run in Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/josephrich98/scrnaseq_packages_and_versioning/blob/main/analysis/ipynb/Seurat_v_Scanpy.ipynb)

Seurat_v_Seurat: [![Run in Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/josephrich98/scrnaseq_packages_and_versioning/blob/main/analysis/ipynb/Seurat_v_Seurat.ipynb)

Seurat_version_comparison: [![Run in Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/josephrich98/scrnaseq_packages_and_versioning/blob/main/analysis/ipynb/Seurat_version_comparison.ipynb)