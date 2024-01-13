Suprachiasmatic Nucleus (SCN) Single-Cell RNA-Sequencing (scRNAseq) Analysis
============================================================================
`KSVelArc`

---

Analysis of SCN scRNAseq data collected from mice. The data were published by [Wen et al. Nat Neurosci 23, 456–467 (2020)](https://doi.org/10.1038/s41593-020-0586-x)

Reference: https://doi.org/10.1038/s41593-020-0586-x \
Data: [GSE117295](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117295), [GSE118403](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118403), [GSE132608](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132608)

<br>

SCRIPTS AND NOTEBOOKS:

`preprocessing.py`
<pre>Python script to process the scRNAseq samples.</pre>

`run_preprocessing.sh`
<pre>Shell script to run <b>preprocessing.py</b></pre>

`model.py`
<pre>Python script to create a scVI model on the data.</pre>

`run_model.sh`
<pre>Shell script to run <b>model.py</b></pre>

`Cluster_ID_Pipeline.ipynb`
<pre>scRNAseq data analysis pipeline for the identification of cell types in the clusters found in the data.</pre>

`scRNAseq.ipynb`
<pre>General information about the scRNAseq technology.</pre>
    
`Wen_etal_2020.ipynb`
<pre>General information about the research published with the data used in this project.</pre>
