# Data pre-processing

This folder contains all the files involving pre-processing and creating of expression dataset into a desired format. It also includes scripts for analysis. 

This desired format is a dictionary of dictionaries. Each dictionary key is a chemical or diseases and the value consist of dictionary with gene and its respective fold changes.

1. [chemicals.py](chemicals.py) notebook consists of pre-processing of chemical datasets.

1. [diseases.py](diseases.py) notebook consists of pre-processing of various disease datasets.

1. [gene_count.py](gene_count.py) notebook for analysis of genes in various datasets

1. [resolver](gene_count.py) notebook for normalizing of differnt namespaces into specific ones (`mondo` for diseases, `pubchem.compound` for chemicals and `ncbigene` for genes)