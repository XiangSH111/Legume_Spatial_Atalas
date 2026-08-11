import scvelo as scv
import loompy
import scvelo as scv
import pandas as pd
import numpy as np
import anndata as ad
import os
import sys
from scipy import sparse as sp

###
adata = scv.read('MX.h5ad', cache=True)
name1 = 'X'
name2 = '../data/MX_metadata'

scv.settings.figdir ='fig_output'
scv.set_figure_params('scvelo')
scv.pp.filter_and_normalize(adata, n_top_genes=3000)
scv.pp.normalize_per_cell(adata)
scv.pp.moments(adata)
scv.pl.proportions(adata,save='_proportion')
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
adata.obsm['X_spatial']=adata.obs[['cx','cy']].values.astype(float)
scv.pl.velocity_embedding_stream(adata, basis='spatial',save='_test')

###
adata.obs=adata.obs.rename(index = lambda x: name1+x)
meta_path = name2
sample_obs = pd.read_csv(os.path.join(meta_path, "cellID_obs.csv"))
cell_umap= pd.read_csv(os.path.join(meta_path, "cell_embeddings.csv"), header=0, names=["Cell ID", "UMAP_1", "UMAP_2"])
cell_clusters = pd.read_csv(os.path.join(meta_path, "cell_clusters.csv"), header=0, names=["Cell ID", "cluster"])

###
sample_one = adata[np.isin(adata.obs.index, sample_obs)]
sample_one_index = pd.DataFrame(sample_one.obs.index)
sample_one_index = sample_one_index.rename(columns = {0:'Cell ID'})
clusters_ordered = sample_one_index.merge(cell_clusters, on = "Cell ID")
umap_ordered = sample_one_index.merge(cell_umap, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
clusters_ordered = clusters_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values
sample_one.obs['clusters'] = clusters_ordered.values

###
cell_celltype = pd.read_csv(os.path.join(meta_path, "cell_celltype.csv"), header=0, names=["Cell ID", "celltype"])
celltype_ordered = sample_one_index.merge(cell_celltype, on = "Cell ID")
celltype_ordered = celltype_ordered.iloc[:,1:]
sample_one.obs['celltype'] = celltype_ordered.values

###
adata = sample_one
###save data
adata.write('Allcelltype_dynamicModel.h5ad', compression = 'gzip')
ident_colours=['#226634','#cdcd9b','#b13824','#e2ba2f','#92b737','#2083aa','#2d488d','#633678']

scv.pl.velocity_embedding_stream(adata, basis='X_umap',size=120,color = "celltype", palette = ident_colours,save="_UMAP",sort_order=True)
scv.pl.velocity_embedding_stream(adata, basis='spatial',color = "celltype", palette = ident_colours,size = 200,alpha =0.7,save="_spatial.pdf",colorbar=True,figsize = (5,8),density = 2.5)
scv.pl.velocity_embedding_grid(adata, basis='spatial', color='celltype', save='spatial_velo.pdf', title='', scale=2,figsize = (5,8),alpha =0.7)


