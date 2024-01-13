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

DIR = sys.argv[1]

adata = sc.read_h5ad(DIR+'.h5ad')

adata.layers['counts'] = adata.X.copy()

# Normalize the counts to 10,000 per cell and log convert the data
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
adata.raw = adata

# Setup the object for scVI analysis, specifying the layer for modeling, categorical covariates, and continuous covariates
scvi.model.SCVI.setup_anndata(adata, layer = 'counts',
                             categorical_covariate_keys=['Sample'],
                             continuous_covariate_keys=['pct_counts_mt', 'total_counts', 'pct_counts_ribo'])

# Initialize an scVI model to capture the variability in the gene expression data
model = scvi.model.SCVI(adata)

start_t = time.time()
# Train the scVI model
model.train() # SLOW
end_t = time.time()

print( (end_t - start_t)/60 )

adata.write_h5ad(DIR+'_integrated.h5ad')
model.save(DIR+'_model')
