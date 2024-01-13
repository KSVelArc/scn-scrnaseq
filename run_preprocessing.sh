#!/bin/bash

ribo_url="http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
dir='data/GSE117295_RAW/'
out='data/wen2020_full.h5ad'


python preprocessing.py "$ribo_url" "$dir" "$out"
