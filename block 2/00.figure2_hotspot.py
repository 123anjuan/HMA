##################################
#hotspot fig2 fiber

##########################
#1,keep protein-coding gene
#2.select top5000 features for co-expression network


import scanpy as sc
import os
import glob
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import pickle
import hotspot
args = sys.argv

#final
adata_3 = sc.read_h5ad(args[1])
adata_3

hs = hotspot.Hotspot(adata_3, model='danb', layer_key = 'counts',latent_obsm_key='X_umap')

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=30,
)
hs_results = hs.compute_autocorrelations(jobs=1)

OUTDIR = "var_5000"
NAME = 'fiber_var5000_nofilter'
HS_RESULTS = ''.join([OUTDIR,"/",NAME,"_hs_results.p"])
LCZ = ''.join([OUTDIR, "/", NAME, "_lcz.p"])
MODULES = ''.join([OUTDIR, "/", NAME, "_modules.p"])
HOTSPOT = ''.join([OUTDIR, "/", NAME, "_hotspot.p"])

with open(HS_RESULTS, "wb") as f:
    pickle.dump(hs_results,f)


#select the genes with significant spatial autocorrelation
hs_genes = hs_results.index[hs_results.FDR < 0.05]

# Compute pair-wise local correlations between these genes
lcz = hs.compute_local_correlations(hs_genes, jobs=40)

with open(LCZ, "wb") as f:
    pickle.dump(lcz,f)
    
###设置module数

#modules = hs.create_modules(
#    min_gene_threshold=10, core_only=True, fdr_threshold=0.05
#)

modules = hs.create_modules(
    min_gene_threshold=10, core_only=False, fdr_threshold=0.05
)
modules.value_counts()

with open(MODULES, "wb") as f:
    pickle.dump(modules, f)

with open(HOTSPOT, "wb") as f:
    pickle.dump(hs,f,protocol=4)    
    
#modules = hs.create_modules(
#    min_gene_threshold=5, core_only=False, fdr_threshold=0.05
#)


results = hs.results.join(hs.modules)

module_scores = hs.calculate_module_scores()

plt.rcParams['figure.figsize'] = (15.0, 12.0)
hs.plot_local_correlations()
