import numpy as np
import pandas as pd
import scanpy as sc
import scvi

import os
import anndata
import matplotlib.pyplot as plt

#Reload harmonized scVI data for placenta:
#results_scvi= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_scVI_harmonization/scvi_integration_1803/SP014_SP082_SP136_villi_scvi_lv15.h5ad"
results_scvi= "/dh-projects/preeclampsia_2022/analysis/spatial_data/SP014_SP082_SP136_placenta_annotations_250322.h5ad"
adata_villi= sc.read_h5ad(results_scvi)
adata_villi.obs['leiden_subclusters_refined02'].cat.categories

#Here, we initialize and fit an AmortizedLDA model on the dataset. 
#Assume that #topic= #cell-types or #leidens (can be tailored according to needs)
#LDA uses the raw count slot of the anndata object (adata.layers)
n_topics = 15 
scvi.model.AmortizedLDA.setup_anndata(adata_villi, layer = "counts")
model = scvi.model.AmortizedLDA(adata_villi, n_topics = n_topics)

#Visualize learned topics:
#By calling model.get_latent_representation(), the model will compute a Monte Carlo estimate of the topic proportions for each cell. 
#Since we use a logistic-Normal distribution to approximate the Dirichlet distribution, the model cannot compute the analytic mean. 
topic_prop = model.get_latent_representation()
topic_prop.head() 

#Save topic proportions in obsm and obs columns.
adata_villi.obsm["X_LDA"] = topic_prop
for i in range(n_topics):
  adata_villi.obs[f"LDA_topic_{i}"] = topic_prop[[f"topic_{i}"]]
  
 
#Color UMAP by topic proportions:
#By coloring by UMAP by topic proportions, we find that the learned topics are generally dominant in cells close together in the UMAP space.
#Sometimes a topic is dominant in multiple clusters in the UMAP, which indicates similarilty between these clusters despite being far apart in the plot.
# Save UMAP to custom .obsm field.
adata_villi.obsm["raw_counts_umap"] = adata_villi.obsm["X_umap"].copy()
sc.pl.embedding(adata_villi, "raw_counts_umap", color = [f"LDA_topic_{i}" for i in range(n_topics)], frameon=False)

#Top features contributing to a topic. 
feature_by_topic = model.get_feature_by_topic()
feature_by_topic.head()

rank_by_topic = pd.DataFrame()

for i in range(n_topics):
    topic_name = f"topic_{i}"
    topic = feature_by_topic[topic_name].sort_values(ascending=False)
    rank_by_topic[topic_name] = topic.index
    rank_by_topic[f"{topic_name}_prop"] = topic.values
	
#Use the same script for decidua but just change the #topics= 25 to model 25 clusters and find topic signatures. 