import os
from src.project_setup import project_setup_function
from src.cellranger_count import cellranger_count_function, cellranger_install_transcriptome_function
from src.read_counts import read_counts_function
from src.kb import kb_ref_function, kb_count_function
from src.downsample_fastqs import downsample_fastqs_function
from scripts.install_cellranger import install_cellranger_function

class FastqProcessor():
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        if not self.data_name:
            self.data_name = self.fastq_link.split("/")[-2]  # 5k_pbmc_v3
        
        self.main_directory = os.path.dirname(os.path.abspath(__file__))
        self.root_directory = os.path.dirname(self.main_directory) # gets the parent directory of the main directory
        self.data_directory = os.path.join(self.main_directory, f"data/{self.data_name}")
        self.original_fastq_directory = os.path.join(self.data_directory, f"{self.data_name}_fastqs")
        self.downsampled_fastq_directory = os.path.join(self.data_directory, "downsampled_fastqs")
        self.output_directory = os.path.join(self.main_directory, f"count_matrix_collection/{self.data_name}")
        self.package_path = os.path.join(self.root_directory, "opt")
        self.reference_path = os.path.join(self.root_directory, "reference_files")
        
    def project_setup(self):
        project_setup_function(self)

    def downsample_fastqs(self):
        downsample_fastqs_function(self)

    def read_counts(self):
        read_counts_function(self)
    
    def kb_ref(self):
        kb_ref_function(self)
        
    def kb_count(self, baseline=False, threads=8):
        kb_count_function(self, baseline=baseline, threads=threads)
        
    def cellranger_install_transcriptome(self):
        cellranger_install_transcriptome_function(self)
    
    def cellranger_count(self, baseline):
        cellranger_count_function(self, baseline=baseline)

    def install_cellranger(self):
        install_cellranger_function(self)
    