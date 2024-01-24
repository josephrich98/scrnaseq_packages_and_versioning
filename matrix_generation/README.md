# to just make conda environment
conda env create -f environment.yml

# change all hyperparameters in config.yaml

# minimal setup (requires fastq files, reference transcriptome (cellranger) or index/t2g (kb) files in reference_files/, and cellranger/kb-python and seqtk installed, all in the structure described below):
python3 main.py

# include all data downloads:
python3 main.py -a

# can also pick and choose which functions to run by running the appropriate file in src

### note: for baseline (no downsampling), notation is seed=0, frac=1

### note: cellranger is only supported on linux AMD architecture, and will not work on the ARM-based docker image

# Expected file structure

~/<root>/  
├── opt/  
│   ├── cellranger-X.X.X/  
│   └── seqtk  
├── reference_files/  
│   ├── cellranger/  
│   │   └── <date_of_install>/  
│   │       └── grch38_transcriptome/  
│   │           └── refdata-gex-GRCh38-2020-A/  
│   │               └── etc.  
│   └── kb/  
│       └── <date_of_install>/
│           ├── index.idx/
│           ├── t2g.txt
│           └── transcriptome.fa
│
└── matrix_generation/
    ├── README.md
    ├── main.py
    ├── conda_environments
    ├── scripts/
    ├── src/
    ├── tests/
    ├── data/
    │   └── <data_name>/
    │       ├── fastqs/
    │       │   ├── file1.fastq.gz
    │       │   ├── file2.fastq.gz
    │       │   └── etc.
    │       └── downsampled_fastqs/
    │           ├── file1.fastq
    │           ├── file2.fastq
    │           └── etc.
    └── count_matrix_collection/
        └── <data_name>/
            └── <matrix_source>/
                ├── <matrix_generation_method>_bc_matrix_collection
                │   └──<FINAL_MATRICES_HERE>
                └── <matrix_generation_method>/
                    └── <additional_parent_file(s)>
                        ├── <matrix>.mtx
                        ├── <genes>.tsv
                        └── <barcodes>.tsv


If preferring to run commands manually:
- build conda environment
conda create -n my_env python=3.11
conda activate my_env
pip install kb-python==0.28.2 gget==0.28.2 pyyaml==6.0.1

- download dataset

- OPTIONAL: downsample fastq files (seqtk)
seqtk sample -s SEED FASTQ_FILE FRACTION | gzip > OUTPUT_FILE_PATH"

- download gget material
gget ref -o ENSEMBL_PATH/info.json -r ENSEMBL_RELEASE -d -w dna,gtf homo_sapiens

- create kb reference transcriptome
kb ref -i index.idx -g t2g.txt -f1 transcriptome.fa FASTA_FILE GTF_FILE

- kb count 
kb count -i index.idx -g t2g.txt -x 10xv3 -o OUTPUT_FOLDER FASTQ1 FASTQ2 ...

- cellranger count
cellranger count --id=SEED_FRAC --transcriptome=TRANSCRIPTOME_PATH --fastqs=FASTQ_DIRECTORY --sample=DATA_NAME --localcores=8 --localmem=64

