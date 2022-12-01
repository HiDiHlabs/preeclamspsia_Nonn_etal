#!/usr/bin/env python
# coding: utf-8
#Author: Olivia Debnath 
#Computational Oncology/BIH (Dr. Naveed Ishaque & Prof. Dr. Roland Eils)

import numpy as np
import pandas as pd
import scanpy as sc

#SP136 cohort samples: processed using 10X V3. 
#Read the Cellbender corrected data for SP136 PE new sample (Donor-13; 073-d):
adata_d13 = sc.read_10x_h5('/data/analysis/preeclampsia_2019/placenta_atlas_2022/SP136_decidua_matrix/P1328_SP136_013/outs/cellbender_filtered.h5')
adata_d13.var_names_make_unique() 

adata_d13.obs['donor_id']= 'Donor13_SP136'
adata_d13.obs['time']= 'late_term'
adata_d13.obs['disease']= 'C'
adata_d13.obs['tissue']= 'decidua'
adata_d13.obs['smoking']= 'No'
adata_d13.obs['cohort']= 'SP136'
adata_d13.obs['library']= '10X v3'
adata_d13.obs['procurement']= 'Oslo'
adata_d13.obs['placental_volume']= 'NA' #rewrite later when info is available. 
adata_d13.obs['gestational_weeks']= 39 
adata_d13.obs['gestational_days']= 277 
adata_d13.obs['maternal_BMI']= 25.013521
adata_d13.obs['maternal_age']= 25

#Read the Cellbender corrected data for SP136 PE new sample (Donor-14; 100-d):
adata_d14 = sc.read_10x_h5('/data/analysis/preeclampsia_2019/placenta_atlas_2022/SP136_decidua_matrix/P1328_SP136_014/outs/cellbender_filtered.h5')
adata_d14.var_names_make_unique() 
adata_d14.obs['donor_id']= 'Donor14_SP136'
adata_d14.obs['time']= 'late_preterm'
adata_d14.obs['disease']= 'PE'
adata_d14.obs['tissue']= 'decidua'
adata_d14.obs['smoking']= 'No'
adata_d14.obs['cohort']= 'SP136'
adata_d14.obs['library']= '10X v3'
adata_d14.obs['procurement']= 'Oslo'
adata_d14.obs['placental_volume']= 'NA'
adata_d14.obs['gestational_days']= 184  
adata_d14.obs['gestational_weeks']= 26 
adata_d14.obs['maternal_BMI']= 20.919421
adata_d14.obs['maternal_age']= 32


#Donor-15 (102-d):
adata_d15 = sc.read_10x_h5('/data/analysis/preeclampsia_2019/placenta_atlas_2022/SP136_decidua_matrix/P1328_SP136_015/outs/cellbender_filtered.h5')
adata_d15.var_names_make_unique() 
adata_d15.obs['donor_id']= 'Donor15_SP136'
adata_d15.obs['time']= 'late_preterm'
adata_d15.obs['disease']= 'PE'
adata_d15.obs['tissue']= 'decidua'
adata_d15.obs['smoking']= 'No'
adata_d15.obs['cohort']= 'SP136'
adata_d15.obs['library']= '10X v3'
adata_d15.obs['procurement']= 'Oslo'
adata_d15.obs['placental_volume']= 'NA'
adata_d14.obs['gestational_days']= 192 
adata_d15.obs['gestational_weeks']= 27 
adata_d15.obs['maternal_BMI']= 27.815882
adata_d15.obs['maternal_age']= 30



#Concatenate the SP136 samples to an anndata:
adata= adata_d13.concatenate(adata_d14, adata_d15)
adata.obs['donor_id'].value_counts()


#Compute QC of the SP136 samples: 
#Select mitochondrial genes: 
adata.var['MT_genes'] = adata.var_names.str.startswith('MT-') 
#ribosomal genes: 
adata.var['Ribo_genes'] = adata.var_names.str.startswith(("RPS","RPL"))
#hemoglobin genes: 
adata.var['HB_genes'] = adata.var_names.str.contains(("^HB[^(P)]"))

#Calculate the QC using scanpy's in built function: 
sc.pp.calculate_qc_metrics(adata, qc_vars=['MT_genes','Ribo_genes','HB_genes'], percent_top=None, inplace=True)
#Percentage of MT transcripts:
mito_genes = adata.var_names.str.startswith('MT-')
#For each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mt2'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
#Add the total counts per cell as observations-annotation to adata:
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
#Percentage of Ribosomal transcripts:
ribo_genes= adata.var_names.str.startswith(("RPS","RPL"))
#For each cell compute fraction of counts in ribo genes vs. all genes
#the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_Ribo2'] = np.sum(adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
#Percentage of HB transcripts:
hb_genes= adata.var_names.str.contains(("^HB[^(P)]"))
#For each cell compute fraction of counts in ribo genes vs. all genes
#the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_HB2'] = np.sum(adata[:, hb_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1


#Before filtering: 
#The #genes per sample look substantially higher than that of 10X V2 samples (SP014)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'donor_id', rotation= 90)


#Filtering: A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. Here we will only consider cells with at least 200 detected genes and genes need to be expressed in at least 3 cells. Please note that those values (QC) are highly dependent on the sample/library preparation method used & also, on the tissue.
#Reference notebook (keep consistency with SP014): SP014_decidua_cellbender_QC_270122.ipynb


adata_new= adata.copy() #make a copy for further filtering.

#Perform all filtering on adata_new: 
sc.pp.filter_cells(adata_new, min_genes=200)
sc.pp.filter_genes(adata_new, min_cells=3)
print(adata_new.n_obs, adata_new.n_vars)

#Fraction of counts assigned to each gene over all cells: 
sc.pl.highest_expr_genes(adata_new, n_top=30) #MALAT1 is an ubiquitous gene. 


#filter for percent mito: keep 5% cut-off. (stick to a strict threshold as we deal with snRNA-seq data)
#The old data will be reanalyzed with 5% cut-off. 
adata_mito_filtered = adata_new[adata_new.obs['pct_counts_MT_genes'] < 5, :]
print(adata_mito_filtered.n_obs) 


#After filtering for the MT genes, check fraction counts: 
sc.pl.highest_expr_genes(adata_mito_filtered, n_top=30) #MALAT1 is an ubiquitous gene. 


#After filtering: 
sc.pl.violin(adata_mito_filtered, ['pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'donor_id', rotation= 90)

#Calculate the XIST counts:
adata_mito_filtered.obs["XIST-counts"] = adata_mito_filtered.X[:,adata_mito_filtered.var_names.str.match('XIST')].toarray()
sc.pl.violin(adata_mito_filtered, ['XIST-counts'],jitter=0, groupby = 'donor_id', rotation= 90)


#Make scatter plots of total counts vs #genes.
#Recommendation from Vienna colleagues: Don't blindly remove cells by "total counts" => might remove the tetraploid decidual EVTs. 
sc.pl.scatter(adata_mito_filtered, x='total_counts', y='pct_counts_MT_genes')
sc.pl.scatter(adata_mito_filtered, x='total_counts', y='n_genes_by_counts')


#Cell-cycle genes: 
#Download from Github: https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
cell_cycle_genes = [x.strip() for x in open('/data/analysis/preeclampsia_2019/herse4olivia/regev_lab_cell_cycle_genes.txt')]
print(len(cell_cycle_genes))

#Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata_mito_filtered.var_names]
print(len(cell_cycle_genes))


adata_final= adata_mito_filtered.copy()
#Save the raw data safely in the 'raw' slot. 
adata_final.raw = adata_final


#Normalize to depth 10,000: 
sc.pp.normalize_per_cell(adata_final, counts_per_cell_after=1e4)
#take the logarithm:
sc.pp.log1p(adata_final)


#Score by cell-cycle: 
sc.tl.score_genes_cell_cycle(adata_final, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pl.violin(adata_final, ['S_score', 'G2M_score'],jitter=0, groupby = 'donor_id', rotation=90)

#Make a count layer to securely store the raw data to proceed with scVI harmonization: 
adata_final.layers["counts"] = adata_final.raw.X.copy()

#Because we still don't the annotations & would like to predict using CellTypist and scANVI transfer annot: 
adata_final.obs['celltypist_majority_voting']= 'Unknown'
adata_final.obs['celltypist_predicted_labels']= 'Unknown'
adata_final.obs['celltypist_majority_voting']



#Save the newly filtered & log-normalized SP136 object for further cell-typing:  
latest= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_anndata/SP136_PE_decidua_cellbender_normalized_matrix_280222.h5ad"
adata_final.write(latest) 


#Read the previously processed SP014 decidua data:
#Each sample is corrected for ambient RNA & random barcode-swapping using Cellbender. 
latest= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_anndata/Decidua_cellbender_normalized_matrix_270122.h5ad"
ldata= sc.read_h5ad(latest)
ldata.obs['disease'].value_counts()
#Concatenate ldata & adata_final:
adata_new= ldata.concatenate(adata_final)
adata_new.obs['donor_id'].value_counts() #after preprocessing. 
#sc.pl.violin(adata_new, ['n_genes_by_counts', 'total_counts', 'pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'disease', rotation= 90)
#sc.pl.violin(adata_new, ['n_genes_by_counts', 'total_counts', 'pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'time', rotation= 90)
#sc.pl.violin(adata_new, ['n_genes_by_counts', 'total_counts', 'pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'procurement', rotation= 90)


#Donor_id(s) later renamed as per manuscript. 
adata_new.obs['donor_id']= adata_new.obs['donor_id'].astype('category')
adata_new.obs['donor_id_reordered']= adata_new.obs['donor_id'].cat.reorder_categories(['Donor-17025-decidua',
       'Donor-SEKR2-decidua', 'Donor-SOZE1-decidua', 'Donor-327-decidua', 'Donor-328-decidua', 'Donor-372-decidua', 'Donor13_SP136', 
        'Donor-274-decidua', 'Donor-389-decidua', 'Donor-419-decidua', 'Donor14_SP136', 'Donor15_SP136'])


plt.rcParams["figure.figsize"] = (7,6)
plt.rcParams["axes.grid"] = False

sc.pl.violin(adata_new, ['n_genes_by_counts'],jitter=0, groupby = 'donor_id_reordered', rotation= 90, save= '_decidua_nGenes_VP.pdf')
sc.pl.violin(adata_new, ['total_counts'],jitter=0, groupby = 'donor_id_reordered', rotation= 90, save= '_decidua_UMI_VP.pdf')
sc.pl.violin(adata_new, ['pct_counts_MT_genes'],jitter=0, groupby = 'donor_id_reordered', rotation= 90, save= '_decidua_pctMT_VP.pdf')
sc.pl.violin(adata_new, ['log1p_n_genes_by_counts'],jitter=0, groupby = 'donor_id_reordered', rotation= 90, save= '_log1p_nGenes.pdf')
sc.pl.violin(adata_new, ['log1p_total_counts'],jitter=0, groupby = 'donor_id_reordered', rotation= 90, save= '_log1p_totalcounts.pdf')

#Scan for NaN or infinity values. 
X = adata_new.X
X = np.log1p(X)
X = np.expm1(X)

mean = X.mean(axis=0).A1
mean[mean == 0] = 1e-12
mean = np.log1p(mean)
np.any(mean == np.inf)

#Read the decidua annotations predicted by Celltypist LR classifier: performed on SP014 (10X V2 were sequenced early)
results_annot= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/placenta_celltypist_model/Decidua_atlas_celltypist_annotated_2701.h5ad" 
adata_annot= sc.read_h5ad(results_annot)

#Transfer the celltypist predicted_labels & majority voting: 
ldata.obs['celltypist_majority_voting']= adata_annot.obs['majority_voting']
ldata.obs['celltypist_predicted_labels']= adata_annot.obs['predicted_labels']

#Repeat the merging: concatenate ldata & adata_final again:
adata_new= ldata.concatenate(adata_final)
adata_new.obs['donor_id'].value_counts() #after preprocessing. 
adata_new.obs['donor_id']= adata_new.obs['donor_id'].astype('category')

adata_new.obs['donor_id_reordered']= adata_new.obs['donor_id'].cat.reorder_categories(['Donor-17025-decidua',
       'Donor-SEKR2-decidua', 'Donor-SOZE1-decidua', 'Donor-327-decidua', 'Donor-328-decidua', 'Donor-372-decidua', 'Donor13_SP136', 
        'Donor-274-decidua', 'Donor-389-decidua', 'Donor-419-decidua', 'Donor14_SP136', 'Donor15_SP136'])

adata_new.obs['donor_id_reordered'].cat.categories

#10X V3 (or, SP136) falls under "Unknown" (needs to be predicted)
#Basis of prediction on 10X V2: Decidua_cellpist_prelims_180222.ipynb 
adata_new.obs['celltypist_predicted_labels'].value_counts() 

adata_new.raw= adata_new #preserve the raw counts safely
#Log-normalize the data: 
sc.pp.normalize_total(adata_new, target_sum=1e4)
sc.pp.log1p(adata_new)

#Subset top 6000 genes (since 4000 HVG couldn't separate the EpCs nicely from NKT, we increased the #HVF).
#All genes worked best: Decidua_celltyping_180222.ipynb
sc.pp.highly_variable_genes(adata_new, n_top_genes= 6000, batch_key= "donor_id", subset=True)

#global seed set to 0 for reproducibility 
import scvi 

#Re-scan for NaN/infinity values:
X = adata_new.X
X = np.log1p(X)
X = np.expm1(X)
mean = X.mean(axis=0).A1
mean[mean == 0] = 1e-12
mean = np.log1p(mean)
np.any(mean == np.inf)

#Input seed labels (for semi-annotated dataset)
#Would be used for refence mapping using scANVI. 
adata_new.obs['seed_labels']= adata_new.obs['celltypist_predicted_labels']

#Transfer of annotation with scANVI: 
#As in the harmonization notebook, we need to register the AnnData object for use in scANVI. Along with batch parameter,  we will give the seed labels for scANVI to use.
scvi.data.setup_anndata(adata_new, layer="counts", batch_key= "donor_id",
                       categorical_covariate_keys= ["time", "gestational_weeks"],
                       continuous_covariate_keys=["n_genes_by_counts", "total_counts", "pct_counts_MT_genes", "XIST-counts"],
                       labels_key="seed_labels")


#Use latent variables=10 & ZINB modelling:
#LV is fin-tuned and checked for n_latent=15, 20, 30 (not shown here)
scvi_model = scvi.model.SCVI(adata_new, n_latent=10, n_layers=2)
scvi_model #details of the pre-initialized model
scvi_model.train(100) #Train the model for 100 epochs. 
scvi_model #trained & converged. 


# Parameters
# scvi_model : SCVI
# Pretrained scvi model
# labels_key : Optional[str] (default: None)
# key in adata.obs for label information. Label categories can not be different if labels_key was used to setup the SCVI model. If None, uses the labels_key used to setup the SCVI model. If that was None, and error is raised.
# unlabeled_category : str
# Value used for unlabeled cells in labels_key used to setup AnnData with scvi.
#Now we can train scANVI and transfer the labels!
#For reference: https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCANVI.from_scvi_model.html
#Here, we need to predict the annotations of "Unknown" (or, unlabeled cells of SP136)
#Belongs to "supervised" analysis. 
scanvi_model = scvi.model.SCANVI.from_scvi_model(scvi_model, 'Unknown')
scanvi_model.train(50)
scanvi_model #scANVI model trained. 

#save the reference model for future use: 
dir_path = "decidua_ref_model/"
scanvi_model.save(dir_path, overwrite=True)

#Now we can predict the missing cell types ("Unknown"), and derive the latent space as normally done in scVI:
adata_new.obs["C_scANVI"] = scanvi_model.predict(adata_new)
#Fetch latent representation of cells: 
adata_new.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata_new)
#Compute KNN graph using "X_scANVI" latent representation:
sc.pp.neighbors(adata_new, use_rep="X_scANVI")
#Compute UMAP using a fixxed random seed (for initialization): 
sc.tl.umap(adata_new, random_state=0, maxiter=500)


#Temporary color-codes. 
adata_new.uns['celltypist_predicted_labels_colors']= ["#9ac6C5", "#ff0000", "#c0c999", "#f87060", "#dec1ff", 
                                                          "#e60000", "#bfff80", "#31cb00", "#fe6776", "#cc33ff", "#ca6680", 
                                                          "#713e5a", "#bf3100", "#d76a03", "#ec9f05", "#8ea604", "#998650", 
                                                          "#e01a4f", "#7a9cc6", "#bde4a7", "#009900"]


#Fine-tune the min_dist=0.3 
sc.tl.umap(adata_new, random_state=0, min_dist= 0.3, maxiter=500)
plt.rcParams["figure.figsize"] = (4,4)
sc.pl.umap(adata_new, color=['seed_labels', 'C_scANVI'])
adata_new.obs['C_scANVI'].value_counts() #distribution per cell-type  


#Read the latest celltypist LR classifier results:
celltye_results= "./decidua_celltypist_new/SP014_SP136_decidua_celltypist.h5ad"
adata_annot= sc.read_h5ad(celltye_results)

#Transfer all celltypist predicted data: 
adata_new.obs['majority_voting'] = adata_annot.obs['majority_voting']
adata_new.obs['predicted_labels']= adata_annot.obs['predicted_labels']
adata_new.obs['conf_score']= adata_annot.obs['conf_score']
adata_new.obs['over_clustering']= adata_annot.obs['over_clustering']

#conf_score: indicates confidence of predicting annotations (1= highest conf & 0= lowest conf)
sc.pl.umap(adata_new, color=['total_counts', 'pct_counts_MT_genes', 'conf_score'], ncols=1)


#Map the legends on data: 
sc.pl.umap(adata_new, color=['C_scANVI', 'predicted_labels', 'majority_voting'], ncols=1, legend_loc= "on data",  
           legend_fontweight= 'normal', legend_fontsize= 'x-small')


#Perform Leiden subclustering: res= 1 (initial)
#Unsupervised analysis: basis of mapping and finalizing clusters (scANVI/celltypist predictions are for validating Leidens)
sc.tl.leiden(adata_new, resolution=1)

sc.pl.umap(adata_new, color=['leiden'], ncols=1, legend_loc= "on data",  
           legend_fontweight= 'normal', legend_fontsize= 'x-small')



#Do a Naive Bayes differential expression analysis:
de_leiden = scanvi_model.differential_expression(groupby= "leiden")
#Save the leiden markers:
de_leiden.to_csv('./decidua_markers_scvi/Decidua_leiden_Bayes_factor_010322.csv')


#Dot-plot visualization of top Leiden res=1 markers: 
markers = {}
cats = adata_new.obs.leiden.cat.categories

for i, c in enumerate(cats):
    cid = "{} vs Rest".format(c)
    cell_type_df = de_leiden.loc[de_leiden.comparison == cid]

    cell_type_df = cell_type_df[cell_type_df.lfc_mean > 1]

    cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3]
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1]

    markers[c] = cell_type_df.index.tolist()[:6]

sc.tl.dendrogram(adata_new, groupby="leiden", use_rep="X_scANVI")
#Top-6 are shown: 
sc.pl.dotplot(adata_new, markers, groupby='leiden', dendrogram=True,
    color_map="Blues", use_raw=True, standard_scale="var", save= '_decidua_leiden_top6_Bayes.pdf')



#Increase the resolution to 2 & over-cluster:
sc.tl.leiden(adata_new, resolution=2, key_added= 'leiden_res2')

sc.pl.umap(adata_new, color=['leiden_res2'], ncols=1, legend_loc= "on data",  
           legend_fontweight= 'normal', legend_fontsize= 'x-small')


#Do a Naive Bayes differential expression analysis again using "leiden_res2" as key:
de_leiden_res2 = scanvi_model.differential_expression(groupby= "leiden_res2")

#Store the top markers found using res=2 labels:  
markers = {}
cats = adata_new.obs.leiden_res2.cat.categories
for i, c in enumerate(cats):
    cid = "{} vs Rest".format(c)
    cell_type_df = de_leiden_res2.loc[de_leiden_res2.comparison == cid]

    cell_type_df = cell_type_df[cell_type_df.lfc_mean > 1]

    cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3]
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1]

    markers[c] = cell_type_df.index.tolist()[:15]

pd.DataFrame(markers).to_csv('Decidua_leiden2_Bayes_top15.csv') #Save as CSV file. 


#Qualitative visualization of decidua integration (scANVI): 
for i in adata_new.obs['time'].cat.categories:
    print(i) 
    fig= sc.pl.umap(adata_new[adata_new.obs['time'] == i], color = 'leiden', return_fig=True, title= i)
    #pdf.savefig(fig)
    

#Split by donor-id
for i in adata_new.obs['donor_id_reordered'].cat.categories:
    print(i) 
    fig= sc.pl.umap(adata_new[adata_new.obs['donor_id_reordered'] == i], color = 'leiden', return_fig=True, title= i)
    #pdf.savefig(fig)

#Separately investigate subclusters found just with 10X V3 samples (for example, two leidens in dMAC compartment)


#Do a Naive Bayes differential expression analysis: on "predicted_labels" (predicted using LR classifier; celltypist)
de_df = scanvi_model.differential_expression(groupby= "predicted_labels")
de_df.to_csv('Decidua_predictedlabels_Bayes_010322.csv')

markers = {}
cats = adata_new.obs.predicted_labels.cat.categories
for i, c in enumerate(cats):
    cid = "{} vs Rest".format(c)
    cell_type_df = de_df.loc[de_df.comparison == cid]

    cell_type_df = cell_type_df[cell_type_df.lfc_mean > 1]

    cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3]
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1]

    markers[c] = cell_type_df.index.tolist()[:15]
    
pd.DataFrame(markers).to_csv('Decidua_predictedlabels_Bayes_top15.csv')

#Try higher cut-off for non-zero proportion of the cell type: 25% or fine-tune parameters as necessary. 

#Barplots for visualizing the Leiden/predicted label composition per donor_id:
#Serves preliminary exploration and not for final analysis in manuscript. 

def make_bar_plots(adata, anno_groupby = 'time', #time: gestational_time (early, late or PE)
                         anno_subsets = 'predicted_labels',
                         anno_pointdef = 'time',
                         save_file = 'Cell_composition_decidua_010322_v1.pdf'):
    
    #Get number of categories: 
    labels = adata.obs[anno_groupby].cat.categories.tolist()
    n_labels = len(labels)
    subset_ids = adata.obs[anno_subsets].cat.categories.tolist()
    n_subsets = len(subset_ids)
    patient_ids = adata.obs[anno_pointdef].unique()
    n_patients = len(patient_ids)

    #Calculate subset fractions per donor: 
    subset_frac = pd.DataFrame(np.empty([n_patients, n_subsets]),
                                    index = patient_ids ,
                                    columns = subset_ids)
    for i in np.arange(n_patients):
        ind2 = adata.obs[anno_pointdef] == patient_ids[i]
        for j in np.arange(n_subsets):
            ind1 = adata.obs[anno_subsets] == subset_ids[j]
            subset_frac.iloc[i,j] = sum(ind1&ind2)
    subset_frac = subset_frac.apply(lambda x: x/sum(x),axis=1)

    #Get donor labels: 
    patient_phenos = pd.DataFrame(index = subset_frac.index, columns=[anno_groupby])
    for i in range(0,len(labels)):
        p_list = adata.obs[anno_pointdef][adata.obs[anno_groupby] == labels[i]].unique().tolist()
        patient_phenos[anno_groupby][p_list] = labels[i]
    patient_phenos['time'] = patient_phenos['time'].astype('category')
    patient_phenos['time'] = patient_phenos['time'].cat.reorder_categories(labels)
    ord_ind = patient_phenos.sort_values('time', ascending=False).index

    subset_frac = subset_frac.loc[ord_ind]

    fig, axes = plt.subplots(1,1, figsize=(8,4))

    subset_frac.plot.barh(stacked=True, grid=False, legend=False, ax=axes, color= ['#023fa5','#7d87b9','#bec1d4','#d6bcc0',
 '#bb7784','#8e063b', '#4a6fe3','#8595e1','#b5bbe3','#e6afb9','#e07b91','#d33f6a','#11c638','#8dd593','#c6dec7',
 '#ead3c6','#f0b98d','#ef9708','#0fcfc0','#9cded6','#d5eae7','#f3e1eb']) #define colors in stacked bar-plots. 

    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    axes.legend(loc='lower right', bbox_to_anchor=(1.5, 0.5), frameon=True)
    sns.despine(left=True, bottom=True, right=True)
    plt.tight_layout()
    plt.savefig(save_file)

import seaborn as sns
make_bar_plots(adata_new, anno_subsets = 'predicted_labels') #celltypist
make_bar_plots(adata_new, anno_subsets = 'leiden_res2') #leidens 

#Subcluster 13, 15 & 18 for demarcating the dMSC vs DSC1 & to separate the dFB1 out. 
sc.tl.leiden(adata_new, resolution=0.35, restrict_to= ('leiden_res2', ['13', '15', '18']), key_added= 'dMSC_leiden')
sc.pl.umap(adata_new, color=['dMSC_leiden'])
adata_new.obs['dMSC_leiden'].value_counts().to_csv('Decidua_leiden_overclust_mesenchymal.csv')

#Do a Bayes differential expression analysis: on "dMSC_leiden" key to find specific subcluster markers. 
de_df = scanvi_model.differential_expression(groupby= "dMSC_leiden")
de_df.to_csv('Decidua_Bayes_dMSC_leiden_010322.csv')
de_df.head()


#Dotplot: top10 markers on the dMSC subcluster labels. 
markers = {}
cats = adata_new.obs.dMSC_leiden.cat.categories

for i, c in enumerate(cats):
    cid = "{} vs Rest".format(c)
    cell_type_df = de_df.loc[de_df.comparison == cid]

    cell_type_df = cell_type_df[cell_type_df.lfc_mean > 1]

    cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3]
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1]

    markers[c] = cell_type_df.index.tolist()[:15]
    
pd.DataFrame(markers).to_csv('Decidua_dMSC_leiden_top15.csv')

sc.tl.dendrogram(adata_new, groupby="dMSC_leiden", use_rep="X_scANVI")
sc.pl.dotplot(adata_new, markers, groupby= 'dMSC_leiden', dendrogram=True,
    color_map="Blues", use_raw=True, standard_scale="var", save= '_decidua_dMSCleiden_top6_Bayes_v1.pdf')


#Predicted labels: how well the dMSC_leiden (mesenchymal subclusters) conform to celltypist predicted labels?! 
df = adata_new.obs.groupby(["predicted_labels", "dMSC_leiden"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
plt.figure(figsize=(12, 8))
_ = plt.pcolor(conf_mat)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted")
plt.ylabel("Observed/Unsupervised")

#Majority voting: do the same confusion matrix using "majority voting" key. 
df = adata_new.obs.groupby(["majority_voting", "dMSC_leiden"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
plt.figure(figsize=(12, 8))
_ = plt.pcolor(conf_mat)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted: Majority vote by LR classifier")
plt.ylabel("Observed/Unsupervised")


#Semi-final merging of cell types/states per Leiden (overclustered):
sc.pl.umap(adata_new, color=['dMSC_leiden'], legend_loc= "on data",  
           legend_fontweight= 'normal', legend_fontsize= 'x-small')

sc.pl.umap(adata_new, color= ['predicted_labels'])


#C38: XCR1, WDFY4, BTLA, SUB1, ZNF366. 
sc.pl.umap(adata_new, color= ['XCR1', 'WDFY4', 'BTLA', 'SUB1', 'ZNF366', 'LGALS2'], ncols=3)

#Maybe, Leiden-24/25 should be subclustered again to better resolve CD14/LYZ monocytes vs dMAC2. 
#Initial mapping of Leidens (based on markers, scANVI & celltypist information)
adata_new.obs['celltype_v1'] = (adata_new.obs["dMSC_leiden"].map(lambda x: {"1": "dNK1", "2": "dNK1", "27": "dNK1", "28": "dNK1",
"34": "dNK_prol", "5": "dNK2", "6": "dTcell", "12": "dTcell", "17": "dTcell", "19": "dTcell", "21": "dTcell", "33": "dTcell", 
"29": "dPlasma_cell", "23": "deported_SCT", "31": "dEVT", "30": "dLEC", "37": "dLEC", "3": "dLEC", "11": "dVEC", "7": "dEpC",
"8": "dEpC", "20": "dEpC", "42": "dEpC", "32": "dEpC", "36": "dGranul", "26": "dMAC1", "16": "dMAC1", "0": "dMAC1", "4": "dMAC1", 
"10": "dMAC1/2 transition", "24": "dMAC2", "25": "dMono_LYZ", "14": "dMAC2", "35": "VWA5A_C35", "22": "dMAC1_sp", "13-15-18,0": "dSMC",
"13-15-18,3": "dSMC/dFB2", "13-15-18,7": "TPM2_dFB2", "13-15-18,2": "DSC2", "13-15-18,6": "DSC2", "13-15-18,1": "DSC1", 
"13-15-18,4": "DSC1", "13-15-18,5": "Uncertain_dMSC", "9": "dFB1", "38": "BTLA_lymphocyte", "39": "C39_VCTs_remove", 
"40": "dFB_donor", "41": "HBG1_C41"}.get(x, x)).astype("category"))

#Visualize the cell annotations so far: 
sc.pl.umap(adata_new, color=['celltype_v1'])


#Filter the C39 VCT contamination plus two other donor specific small clusters (such as L41):
ldata_dec= adata_new[~adata_new.obs['celltype_v1'].isin(['C39_VCTs_remove', 'HBG1_41', 'dFB_donor'])]
sc.pl.umap(ldata_dec, color=['celltype_v1']) #replot the UMAP. 


#Plot dMSC markers for first impression:
#ldata_dec: newly filtered dataset after spurious cluster removal.  
sc.pl.umap(ldata_dec, color=['ZEB1', 'MEIS1', 'SOX5', 'FOXO1', 'GPC6', 'TWIST1', 'TWIST2'])

#Do a Bayes differential expression analysis: "celltype_v1" (refined cluster key)
de_df = scanvi_model.differential_expression(groupby= "celltype_v1")
de_df.to_csv('Decidua_Bayes_celltypev1_010322.csv')


#Plot top6 markers for the predicted_labels: 
markers = {}
cats = ldata_dec.obs.celltype_v1.cat.categories

for i, c in enumerate(cats):
    cid = "{} vs Rest".format(c)
    cell_type_df = de_df.loc[de_df.comparison == cid]

    cell_type_df = cell_type_df[cell_type_df.lfc_mean > 1]

    cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3]
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1]

    markers[c] = cell_type_df.index.tolist()[:6]
    
sc.tl.dendrogram(ldata_dec, groupby= "celltype_v1", use_rep="X_scANVI")

sc.pl.dotplot(ldata_dec, markers, groupby= 'celltype_v1', dendrogram=True,
    color_map="Blues", use_raw=True, standard_scale="var", save= '_decidua_celltypev1_top6_Bayes_v1.pdf')



ldata_new= adata_new.copy()
#ldata_new.obs= ldata_new.obs.astype('category')
#ldata_new.obs['maternal_BMI']= ldata_new.obs['maternal_BMI'].astype('category')

#Convert to strings to save the anndata as H5 class below: 
ldata_new.obs= ldata_new.obs.astype('string')
ldata_new.obs.columns.values.dtype

#Now, h5 class can convert it to categoricals while writing anndata:  
#Note that, new path on Eils-HPC is: "/dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/decidua_anndata/SP014_SP136_decidua_analysis_040322.h5ad"
dec_results= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_anndata/SP014_SP136_decidua_analysis_010322.h5ad"
ldata_new.write(dec_results)

#Read the data above & check:
adata_test= sc.read_h5ad(dec_results)
adata_test.raw.X #dim looks correct. 
sc.pl.umap(adata_test, color=['celltype_v1']) #works! 

del adata_new 

#Part-II: Semi-finalizing the cell annotations for manuscript. Please adhere to steps below. 
#For more detailed explanation on the underlying biology of clusters, have a look at: Decidua_immune_stromal_analysis_exploratory_01032022.ipynb
#Note that, new path on Eils-HPC is: "/dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/decidua_anndata/SP014_SP136_decidua_analysis_040322.h5ad"
#Section below is exerc
dec_results= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_anndata/SP014_SP136_decidua_analysis_010322.h5ad"
ldata_dec01= sc.read_h5ad(dec_results)

#Redo the filtering step: 
ldata_new= ldata_dec01[~ldata_dec01.obs['celltype_v2'].isin(['C39_VCTs_remove', 'HBG1_C41', 'dFB_donor'])] 


plt.rcParams["figure.figsize"] = (4,4)
plt.rcParams["axes.grid"] = False

#Do the renaming of MAC subclusters:
ldata_new.obs['celltype_v3'] = (ldata_new.obs["celltype_v2"].map(lambda x: {"dMAC1/2 transition": "dMAC2", 
"dMAC_sp": "NKT", "dMAC2": "IL17RA_dMAC3", "BTLA_lymphocyte": "dDC",
"dGranul": "dGranulocyte"}.get(x, x)).astype("category"))

sc.pl.umap(ldata_new, color=['celltype_v3'], ncols=1)


#Notes: 
#dMAC1: LYVE1/F13A1 positive mac; also, CD163 high & express other canonical mac-genes such as RBPJ, MRC1, SPP1, VSIG4 & SELENOP. CD209, CD86 & MAF are also expressed.
#dMAC-1/2 transition (can be named as dMAC2): Relatively lower exp of CD163/F13A1 etc; presence of few robust genes like FTH1, GLUL, CTSL & CTSB. The CTS (cathepsin) genes are endopeptidases & are described before in the context of tumor-associated macrophages. They can cleave intercellular adhesion molecules as E-cadherin and JAMs & hence, can enhance migration & tumor progression. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4180975/ & https://www.nature.com/articles/nrc2816
#dMAC2 (can be renamed as dMAC3): Robust for S100A9, LYN, MNDA, IL17RA, JAML, MXD1, S100A8, LUCAT2 & CSF3R. Express TLR2 & DOCK4 but there is a negligible exp of well-known genes like F13A1, MRC1, SPP1, etc. RBPJ is expressed in some cells.

#- CSF3R can provide some clues known to regulate granulopoiesis, neutrophil function, and 
#hematopoietic stem cell mobilization: https://rupress.org/jem/article/208/2/251/40868/Expressionof-the-G-CSF-receptor-in-monocytic

 #- IL17RA is also important & it’s found that IL17 receptors could be upregulated on macrophages in 
#vitro, as well as in vivo, in response to different inflammatory stimuli: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4292791/

#dMAC_sp: Seems like NKT cells. Express NCAM1, CD96, GNLY & then, THEMIS. Around 20% of cells also show NKG7. 
#They also show some macrophage gene markers (not robust though)- so, wondering if there are some contaminations. Later removed because it's ambiguous. 

#Re-label the MAC subclusters a/c to the description above:
ldata_new.obs['celltype_v3'] = (ldata_new.obs["celltype_v2"].map(lambda x: {"dMAC1/2 transition": "dMAC2", 
"dMAC_sp": "NKT", "dMAC2": "IL17RA_dMAC3", "BTLA_lymphocyte": "dDC",
"dGranul": "dGranulocyte"}.get(x, x)).astype("category"))


#Subset the dLEC based leidens:
ldata_lec= ldata_new[ldata_new.obs['dMSC_leiden'].isin(['3', '30', '37'])]

#Plot the main genes from Teichmann's figure-3:
lec_genes= ['STMN1', 'TFF3', 'IFITM3', 'SRP14', 'RAMP2', 'CCL21', 'VAMP5', 'LYVE1', 'MARCKSL1',
           'EFNA5', 'ANO2', 'FABP5']
		   
sc.pl.dotplot(ldata_lec, lec_genes, groupby= 'dMSC_leiden', dendrogram=False,
    color_map="Blues", use_raw=True, standard_scale="var")	

#Do the renaming of LEC subclusters & mark the dLECP (proliferating):
ldata_lec.obs['celltype_v3'] = (ldata_lec.obs["dMSC_leiden"].map(lambda x: {"30": "dLECp", 
"3": "dLEC", "37": "dLEC"}.get(x, x)).astype("category"))
print(ldata_lec.obs['celltype_v3'].value_counts())	#print the composition 

#Fetch the barcodes for LEC/LECP & map back to our main dataset:
lec_cluster = pd.Series(ldata_lec.obs['celltype_v3'])
ldata_new.obs["celltype_v3"]= ldata_new.obs["celltype_v2"] #create a copy. 
ldata_new.obs["celltype_v3"] = ldata_new.obs["celltype_v3"].astype(str)
ldata_new.obs["celltype_v3"].loc[lec_cluster.index] = lec_cluster
ldata_new.obs["celltype_v3"].value_counts() #532 LECP labeled. 

#Repeat the renaming of MAC subclusters: 
ldata_new.obs['celltype_v3'] = (ldata_new.obs["celltype_v2"].map(lambda x: {"dMAC1/2 transition": "dMAC2", 
"dMAC1_sp": "NKT", "dMAC2": "IL17RA_dMAC3", "BTLA_lymphocyte": "dDC",
"dGranul": "dGranulocyte"}.get(x, x)).astype("category"))

#Repeat: 
ldata_new.obs["celltype_v3"] = ldata_new.obs["celltype_v3"].astype(str)

#Fetch the index of the LEC cells: 
ldata_new.obs["celltype_v3"].loc[lec_cluster.index] = lec_cluster


#Subset the stromal compartment (later repeat without dMSC & DSC3):
#Current, dFB2/dFB1 can be improved using the "predicted_labels" annotations from CellTypist:
ldata_stromal= ldata_new[ldata_new.obs['celltype_v3'].isin(['DSC1', 'DSC2', 'dFB1', 'dFB2', 'dSMC', 'dMSC', 'DSC3'])]
print(ldata_stromal.obs['celltype_v3'].value_counts()) 

#Compare both labels: "celltype_v3" is Leiden derived. 
sc.pl.umap(ldata_stromal, color= ['celltype_v3', 'predicted_labels'], ncols=1) 

#Plot dEVT genes: 
sc.pl.umap(ldata_stromal, color= ['MYCNUT', 'NOTUM', 'HPGD', 'LAIR2', 'HLA-G'], ncols=5)

#Since DSC3 cluster express very mixed markers: NOTUM, HPGD, HLA-G along with stromal & macrophage markers, 
#it's very likely biological doublets (contaminations from sampling):
#We've also discussed this separately in: Decidua_immune_stromal_analysis_exploratory_01032022.ipynb 
ldata_new= ldata_new[ldata_new.obs['celltype_v3']!= 'DSC3']
ldata_new.obs['celltype_v3'].cat.categories

#Re-subset stromal cells (after removing DSC3):
ldata_stromal= ldata_new[ldata_new.obs['celltype_v3'].isin(['DSC1', 'DSC2', 'dFB2', 'dSMC', 'dMSC'])]
print(ldata_stromal.obs['celltype_v3'].value_counts()) 

#Look if any dEVT cells are wrongly classified as DSC2:
print(pd.crosstab(ldata_stromal.obs['celltype_v3'], ldata_stromal.obs['predicted_labels'])) 

#Also, this:
print(ldata_stromal.obs['predicted_labels'].value_counts()) 

#Since DSC2 express certain immune genes (agreement with Teichmann's dataset), a part of them were misclassified as immune cells using supervised LR classifiers. 
#Correct/re-map the stromal cell types/states: 
ldata_stromal.obs['predicted_labels_new'] = (ldata_stromal.obs["predicted_labels"].map(lambda x: {"dLEC": "DSC_2",
"dLEC_dysfunctional": "DSC_2", "dMAC_activated": "DSC_2", "dMAC_classical": "DSC_2", "dMonocyte": "DSC_2",   
"dNK_1": "dSMC", "dNK_2": "dSMC", "dNK_prol": "dSMC",  "dTcell": "DSC_2",  "dSCT": "DSC_1",  "dGranulocyte": "DSC_2",
"dEpC": "DSC_2", "dVEC": "dSMC"}.get(x, x)).astype("category"))

#Pedantic:
ldata_stromal.obs['predicted_labels_new']= (ldata_stromal.obs['predicted_labels_new'].map(lambda x: {"DSC_2": "DSC2", 
"DSC_1": "DSC1", "dFB_2": "dFB2", "dFB_1": "dFB1"}.get(x, x)).astype("category"))
print(ldata_stromal.obs['predicted_labels_new'].value_counts()) 


#Finally, store the refined stromal subgroup labels to the parent anndataset i.e., ldata_new: 
stromal_cluster = pd.Series(ldata_stromal.obs['predicted_labels_new'])
ldata_new.obs["celltype_v4"]= ldata_new.obs["celltype_v3"] #create a copy. 
ldata_new.obs["celltype_v4"] = ldata_new.obs["celltype_v4"].astype(str)
#Fetch the index of the stromal cells: 
ldata_new.obs["celltype_v4"].loc[stromal_cluster.index] = stromal_cluster
ldata_new.obs["celltype_v4"].value_counts()

#Visualize the semi-final annotations until now:
sc.pl.umap(ldata_new, color= ['celltype_v4']) 

#Validate the stromal subgroups (latest) with approprite markers:
#Reference notebook: Decidua_immune_stromal_analysis_exploratory_01032022.ipynb 
#Replot with the fig3a_genes (from Teichmann Nature): 
stromal_genes= ['IGF1', 'MME', 'SOX5', 'ZEB1', 'TWIST2', 'SDK1', 'MEIS1', 'SYNPO2', 'ABCG2', 'NANOG',
               'ABCC9', 'WNT5A', 'MCAM', 'PDGFRB','ABI3BP', 'IL15','FOXO1', 'SCARA5', 'PLAGL1', 'ING1', 
                'PRL', 'IGFBP1', 'IGFBP2', 'IGFBP6', 'CXCL14', 'GLRX', 'MGP', 
                'CPEB4', 'ABTB2', 'RGCC', 'SIPA1L2', 'SYT1',  
               'APOD', 'APOE', 'CFD', 'CLU', 'LGALS1', 'LUM', 'S100A6', 'PRL', 'TIMP1', 'TIMP2', 'VIM', 
                'B2M', 'RBP1', 'ACTB', 'TIMP3', 'EPYC','CSKMT', 'IFITM2','COL1A1', 'COL1A2', 'COL3A1', 
                'EDIL3', 'PDLIM3', 'PDLIM5', 'PDLIM7', 'FN1', 'TAGLN', 'SULF1', 'THBS1', 'THBS2', 'TNC', 'TPM2', 
                'SASH1', 'SORBS2', 'SLIT2', 'NOTCH3', 'SPARC', 'CARMN', 'GUCY1A2', 'ACTA2', 'ADGRL3', 'RGS5','EGFLAM', 'MYO1B', 'GUCY1A1', 
                'ADARB2', 'LURAP1L', 'COL18A1', 'CPM', 'CCDC102B', 'EPS8', 'ID2']

#Dot-plot depicting stromal features: 
sc.tl.dendrogram(ldata_stromal, groupby= 'predicted_labels_new') 
sc.pl.dotplot(ldata_stromal, stromal_genes, groupby= 'predicted_labels_new', dendrogram=True, color_map="Blues", use_raw=True, standard_scale="var") 

#Scan for MAC (immune genes) in the stromal subcluster:
bayes2_genes= ['CD163', 'F13A1', 'RBPJ', 'SELENOP', 'SLC9A9', 'C1QB', 'C1QA', 
               'IGFBP1', 'IGFBP2', 'IGFBP4', 'PRL', 'DCN','RGCC', 'TAGLN', 'SPARC',
             'COL1A1', 'COL1A2', 'COL3A1', 'NOTUM', 'MYCNUT', 'HLA-G', 'LAIR2', 'HPGD']

#Plot without the dEVT cluster:
sc.pl.dotplot(ldata_stromal[ldata_stromal.obs['predicted_labels_new']!= 'dEVT'], stromal_genes, groupby= 'predicted_labels_new',  dendrogram=False, color_map="Blues", use_raw=True, standard_scale="var") 

#Subset NKC groups:
ldata_nkc= ldata_new[ldata_new.obs["celltype_v4"].isin(['dNK1', 'dNK2', 'dNK_prol'])]

#Plot few cycling/prol genes to double check dNKprol: 
prol_genes= ['CD96', 'GNLY',  'NCAM1', 'NKG7', 'KLRC1', 'KLRC2', 
             'TOP2A', 'MKI67', 'TPX2', 'MELK', 'CENPF', 'NCAPG']

sc.pl.dotplot(ldata_nkc, prol_genes, groupby= 'dMSC_leiden', dendrogram=False, color_map="Blues", use_raw=True, standard_scale="var")

#Remap categories for dGranulocyte a/c markers:
ldata_new.obs['celltype_v4'] = (ldata_new.obs["celltype_v4"].map(lambda x: {"VWA5A_C35": "dGranulocyte", "dGranulocyte": "dGranul_new"}.get(x, x)).astype("category"))

#Rewrite to a new anndata: 
results= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_anndata/SP014_SP136_decidua_analysis_040322.h5ad"
ldata_new.write(results)

#Subset L24 for monocyte: 
ldata_mono= ldata_new[ldata_new.obs['dMSC_leiden']== '24']
ldata_mono.obs['celltype_v4'] #should be monocytes: CD14, FCN1, LYST etc. 

#Relabel IL17RA_dMAC3: dMono_LYZ (in manuscript, it's dMono1). 
ldata_mono.obs['celltype_v4'] = (ldata_mono.obs['celltype_v4'].map(lambda x: {"IL17RA_dMAC3": "dMono_LYZ"}.get(x, x)).astype("category"))

#Finally, do this to refine the monocyte annotations in "ldata_new": 
mono24_cluster = pd.Series(ldata_mono.obs['celltype_v4'])
ldata_new.obs["celltype_v5"]= ldata_new.obs["celltype_v4"] #create a copy. 
ldata_new.obs["celltype_v5"] = ldata_new.obs["celltype_v5"].astype(str) #set as string. 
#Fetch the index of the monocyte (L24) cells: 
ldata_new.obs["celltype_v5"].loc[mono24_cluster.index] = mono24_cluster
print(ldata_new.obs["celltype_v5"].value_counts()) #almost final annotations. 

#Rewrite the existing anndata to save back the monocyte labels: 
#Note that, current directory on Eils-HPC: /dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/decidua_anndata/SP014_SP136_decidua_analysis_040322.h5ad
results= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_anndata/SP014_SP136_decidua_analysis_040322.h5ad"
ldata_new.write(results)


#Store the latest color codes for category label:
ldata_new.uns['celltype_v5_colors']= ['#33ccff', '#cc3300', '#669999', '#334d4d', '#003366', 
    '#330000','#cc0066', '#e60099', '#ffcc00', '#f2e6ff', '#c65353', '#ff66b3', '#660033',
    '#b3b300', '#4d4d00', '#664400', '#004d00', '#006699', '#80d4ff', '#001a13', '#990000', 
    '#000080', '#3d5c5c', '#9966ff', '#e6005c']
	
#Visualize the UMAP embedding using the latest annotations: 
sc.pl.umap(ldata_new, color= ['celltype_v5'])

#Switch to "/Decidua_celltyping/Decidua_analysis_manuscript_2022.ipynb" to see further steps & generation of decidua related manuscript figures. 




	
			  