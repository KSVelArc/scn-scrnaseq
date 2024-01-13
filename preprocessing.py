#!/usr/bin/env python

import os, sys
import time
import scanpy as sc
import scvi
scvi.settings.seed = 0
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.sparse import csr_matrix
#conda install -c conda-forge numba


def pp(csv_path):
    # Read sample file and convert the dense matrix to a condensed sparse matrix
    adata = sc.read_csv(csv_path).T
    adata.obs_names_make_unique
    adata.X = csr_matrix(adata.X)
    
    # Keep genes expressed in at least 10 cells
    sc.pp.filter_genes(adata, min_cells=10)
    
    # Keep the 2000 most variable genes that best describe the data. subset: , flavor:
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, flavor='seurat_v3')
    
    # DOUBLET REMOVAL
    
    # Train a scVI (single-cell Variational Inference [Lopez et al., 2018]) model. Set up the AnnData object and train the model.
    scvi.model.SCVI.setup_anndata(adata)
    vae = scvi.model.SCVI(adata)
    vae.train()
    
    # Train a Solo (Bernstein et al., 2020) model for doublet detection, passing the pre-trained scVI object. 
    solo = scvi.external.SOLO.from_scvi_model(vae)
    solo.train()
    df = solo.predict()
    df['prediction'] = solo.predict(soft = False)
    
    # Remove the '-0' from the cell barcodes
    #df.index = df.index.map(lambda x: x[:-2])
    
    # Add a column representing the difference between the doublet and singlet predictions
    df['dif'] = df.doublet - df.singlet
    
    # Create a new object containing doublets 
    doublets = df[(df.prediction == 'doublet') & (df.dif > 1)]
    
    # Reload the sample data 
    adata = sc.read_csv(csv_path).T
    adata.obs['Sample'] = csv_path.split('_')[2] #'raw_counts/GSM5226574_C51ctr_raw_counts.csv'
    
    # Add a new boolean column indicating doublet or not
    adata.obs['doublet'] = adata.obs.index.isin(doublets.index)
    adata = adata[~adata.obs.doublet]
    
    # FILTER CELLS
    
    # Remove cells with fewer than 200 genes
    sc.pp.filter_cells(adata, min_genes=200) 
    # Remove genes expressed in fewer than 3 cells
    #sc.pp.filter_genes(adata, min_cells=3)
    
    # REMOVE MITOCHONDRIAL AND RIBOSOMAL GENES
    
    # Create membership columns for mitochondrial and ribosomal genes
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)
    
    # Remove cells containing a threshold of outlier cells, mt and ribo genes
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, 0.98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata = adata[adata.obs.pct_counts_mt < 20]
    adata = adata[adata.obs.pct_counts_ribo < 2]
    
    return adata


# Passed variables
ribo_url = sys.argv[1]
DIR = sys.argv[2]
OUT = sys.argv[3]

ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)


start_time = time.time()


out = []
patterns = ['flash.csv.gz','_CT22.csv.gz']    # 
for file in os.listdir(DIR):
	if any(file.endswith(pattern) for pattern in patterns):
		continue
	else:	
		out.append(pp(DIR + file))

# Concatenate all the objects
adata = sc.concat(out)
    
end_time = time.time()

print( (end_time - start_time)/60 )

# Save to Scanpy format
adata.write_h5ad(OUT)
df.to_csv(OUT+'.tsv', index=False, sep='\t')

# UserWarning: Observation names are not unique. To make them unique, call `.obs_names_make_unique`.
# utils.warn_names_duplicates("obs")
