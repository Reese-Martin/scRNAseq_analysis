# this script is designed to outline a basic qc workflow for scRNAseq data. I will use publicly available date
# for this work, but the major components are interoperable between data sets

import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation

### basic settinngs
sc.settings.verbosity = 0 #only show errors

sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
) # configure figure plotting

# load bone marrow dataset from the 2021 NeurIps conference
adata = sc.read_10x_h5(
    filename="filtered_feature_bc_matrix.h5",
    backup_url="https://figshare.com/ndownloader/files/39546196",
)
adata.var_names_make_unique()
adata
# should see anndata object with the shape n_obs x n_vars = 16934 x 36601 which corresponds to
# barcodes x transcripts (barcodes being unique cells and transcripts being the genes/gene things identified)

### now that data is loaded, we will cull cells that are low quality, cells duos, and cell-free rna contaminated samples
# low quality cells will be selected based on three (mutual) criteria: low counts/barcode, low genes/barcode, MT gene fraction

# first mitochandrial, ribosomal, and heme genes.
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))

# using our added Var entries, now calculate cell qc metrics. Note, this function will add several entries to the
# var and obs sections of our anndata object
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
)

### what are the new fields? some important ones are below
# n_genes_by_counts in obs is gene counts per cell
# total_counts in obs is the total number of counts for a cell. aka library size
# pct_counts_mt is the % of counts that are MT genes. a quick check shows that highest % MT in our data
max(adata.obs["pct_counts_mt"])

# plot time!
p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
p2 = sc.pl.violin(adata, "pct_counts_mt")
p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
