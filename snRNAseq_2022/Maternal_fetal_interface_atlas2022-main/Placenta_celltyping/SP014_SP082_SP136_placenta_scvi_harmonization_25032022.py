#!/usr/bin/env python
# coding: utf-8

#Analysis performed by: Olivia Debnath
#Computational Oncology; Supervisors: Dr. Naveed Ishaque & Prof.Dr.Roland Eils. 
import numpy as np
import pandas as pd
import scanpy as sc

#Set the figure resolution. 
sc.settings.set_figure_params(dpi=80)

#Read the anndata H5AD having the controls plus newly added SP136 PE diseased samples. 
#Important note (after September'22): Old path /data/analysis/preeclampsia_2019/.. is changed to /dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/villi_anndata for storage in Eils-HPC cluster. 
results_placenta= '/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_anndata/SP014_SP082_SP136_placenta.h5ad'
adata= sc.read_h5ad(results_placenta)
adata.layers["counts"]= adata.X #Store the raw data safely to the "count" layer. 
adata.raw= adata #Store the raw data safely.
sc.pp.normalize_total(adata, target_sum=1e4) #Normalize to depth 10,000: 
sc.pp.log1p(adata) #Log-transform the data

adata.raw.X #Check matrix dimension 

#Create a cosmetic copy of the anndata for HVG subselection & data harmonization using scVI. 
adata_hvg= adata.copy()

#Selected first 6000 highly variable features (using donor_ID as a batch key):
#Setting batch key performs light-weight batch correction to prevent selection of very sample specific genes. 
sc.pp.highly_variable_genes(adata_hvg, n_top_genes= 6000, batch_key= "donor_id", subset=True)
adata_hvg.var.index

#Import scvi for data harmonization; global seed is automatically set to 0 for reproducibility. 
import scvi

#Registers anndata object for scVI models.
#A mapping will be created between data fields used by scvi to their respective locations in adata. 
#This step will also compute the log mean and log variance per batch for the library size prior. 
#batch: metadata consisting of donor_id; 
#categorical_covariate_keys: keys in adata.obs that correspond to categorical data. 
#continuous_covariate_keys: keys in adata.obs that correspond to continuous data.  

scvi.data.setup_anndata(adata_hvg, layer= "counts", batch_key= "donor_id",
                       categorical_covariate_keys= ["time", "gestational_days", "library", "procurement"],
                       continuous_covariate_keys=["n_genes_by_counts", "total_counts", "pct_counts_MT_genes", "XIST-counts"])

#Train the variational autoencoder using LV=15 after encoding the covariates: 
#Note that the latent variables are not naturally ordered like the PC(s). 
#Number of nodes per hidden layer= 128 (scVI default)
#n_layers: Number of hidden layers used for encoder and decoder NNs.
#Dropout rate for neural networks= 0.1 (scVI default)
#Dispersion parameter of NB is constant per gene across cell (here, nuclei)
#Zero-inflated negative binomial distribution is used for modeling gene likelihood distribution. 
vae_lv15 = scvi.model.SCVI(adata_hvg, n_layers=2, n_latent=15, gene_likelihood= 'zinb')
vae_lv15.train()
vae_lv15.save("placenta_scvi_lv15/") #Save the trained VAE for future reference.  

#After the training is completed, we can evaluate the latent representation of each cell in the dataset and add it to adata. 
adata_hvg.obsm["X_scVI"] = vae_lv15.get_latent_representation()
adata_hvg.obsm["X_normalized_scVI"] = vae_lv15.get_normalized_expression()

#Finally, we can compute KNN graph and visualize it with UMAP:
sc.pp.neighbors(adata_hvg, use_rep="X_scVI")
sc.tl.leiden(adata_hvg, resolution=2) #over-cluster using Leiden resolution=2. 
sc.tl.umap(adata_hvg, random_state= 0, maxiter=500)

#Visualize the UMAP using following fields: 
sc.pl.umap(adata_hvg, color= ['celltype_v1', 'leiden', 'library', 'time'], ncols=1)

#Qualitatively visualize the integration. Split by gestational time: early, late term controls & preterm diseased. 
for i in adata_hvg.obs['time'].cat.categories:
    print(i) 
    fig= sc.pl.umap(adata_hvg[adata_hvg.obs['time'] == i], color = 'celltype_v1', return_fig=True, title= i+' samples')
    #pdf.savefig(fig)

#Create a cosmetic copy of adata_hvg to check robustness of UMAP (change random_state to 42)
ldata_umap= adata_hvg.copy() 
sc.tl.umap(ldata_umap, random_state=42, maxiter=500) 
sc.pl.umap(ldata_umap, color= ['celltype_v1']) #almost same data structure is obtained. 


#Create a copy again & increase the min_dist to 0.6 to evaluate vSTBjuv proximity to vVECs. 
#Maxiter=1000 for better convergence of UMAP. 
ldata_umap= adata_hvg.copy() 
sc.tl.umap(ldata_umap, random_state=0, maxiter=1000, min_dist=0.6) 
sc.pl.umap(ldata_umap, color= ['celltype_v1'])


#Split the UMAP graph by library chemistry used to prepare samples i.e., 10X V3 (all new samples) vs 10X V2 (old samples aka SP014)
for i in adata_hvg.obs['library'].cat.categories:
    print(i) 
    fig= sc.pl.umap(adata_hvg[adata_hvg.obs['library'] == i], color = 'celltype_v1', return_fig=True, title= i+' samples')
    #pdf.savefig(fig)

#Assign color codes a/c 2021 submission. 
adata_hvg.uns['celltype_v1_colors']= ['#ff0000', '#c0c999', '#9e1a1a', '#721313', '#5c7aff', '#dec1ff', '#8c29ff',
                                 '#fe6776', '#63264a', '#ff99cc', '#ff0080', '#cc33ff', '#a799b7', '#00b3b3', '#004d4d', 
                                 '#fd96A9']


#Reorder the donor_id so that samples are grouped by early, late term controls followed by preterm PE. 
#The donor_id(s) are later renamed a/c to Graz clinical patient_Id (check Placenta_anndata_zenodo_110722.ipynb for reference)
adata_hvg.obs['donor_id_reordered']= adata_hvg.obs['donor_id'].cat.reorder_categories(['Donor-17025-villi', 'Donor-SEKR2-villi', 'Donor-SOZE1-villi',
'Donor11', 'Donor12', 'Donor34', 'Donor37', 'Donor38', 'Donor39', 'Donor41', 'Donor-327-villi', 'Donor-328-villi', 'Donor-372-villi',
'Donor-Term01', 'Donor-Term02', 'Donor-Term03', 'Donor-389-villi', 'Donor-419-villi', 'Donor-557_1-villi', 'Donor-557_2-villi', 'SP136_022v',
'SP136_023v'])

#Rename a/c to the Extended table of the manuscript: 
adata_hvg.obs['donor_id_renamed'] = (adata_hvg.obs["donor_id_reordered"].map(lambda x: {"Donor-17025-villi": "17-025-v",
"Donor-SEKR2-villi": "17-022-v", "Donor-SOZE1-villi": "17-021-v", "Donor11": "18-082-01-v", "Donor12": "18-108-02-v", 
"Donor34": "18-032-v", "Donor37": "18-017-v", "Donor38": "18-033-v", "Donor39": "18-098-v", "Donor41": "20-027-v", 
"Donor-327-villi": "327-v", "Donor-328-villi": "328-v", "Donor-372-villi": "372-v", "Donor-Term01": "TRM1-v", 
"Donor-Term02": "TRM2-v", "Donor-Term03": "TRM3-v", "Donor-389-villi": "389-v", "Donor-419-villi": "419-v", 
"Donor-557_1-villi": "577-1-v", "Donor-557_2-villi": "577-2-v", 
"SP136_022v": "M181-v", "SP136_023v": "PLA120-v"}.get(x, x)).astype("category"))

#Reorder the renamed categories to follow the groupings within category:
adata_hvg.obs['donor_id_reordered']= adata_hvg.obs['donor_id_renamed'].cat.reorder_categories(["17-025-v",
"17-022-v", "17-021-v", "18-082-01-v", "18-108-02-v", "18-032-v", "18-017-v", "18-033-v", "18-098-v", 
"20-027-v", "327-v",  "328-v", "372-v", "TRM1-v", "TRM2-v", "TRM3-v", "389-v", "419-v", 
"577-1-v", "577-2-v", "M181-v", "PLA120-v"])


#Barplots for each patient to visualize the cell type composition:
#Split sample composition by batch= donor_id:
#Reference notebook: Placenta_scVI_controls_exploratory_celltyping_LV15_S3.ipynb
#Visualization inspired by Broad Institute's sepsis project. 
def make_bar_plots(adata, anno_groupby = 'time',
                         anno_subsets = 'celltype_v1',
                         anno_pointdef = 'donor_id_reordered',
                         save_file = 'SP14_SP82_SP136_placenta_atlas_comp_220322.pdf'):
    
    #Get number of categories
    labels = adata.obs[anno_groupby].cat.categories.tolist()
    n_labels = len(labels)
    subset_ids = adata.obs[anno_subsets].cat.categories.tolist()
    n_subsets = len(subset_ids)
    patient_ids = adata.obs[anno_pointdef].unique()
    n_patients = len(patient_ids)

    #Calculate subset fractions for each patient: 
    subset_frac = pd.DataFrame(np.empty([n_patients, n_subsets]),
                                    index = patient_ids ,
                                    columns = subset_ids)
    for i in np.arange(n_patients):
        ind2 = adata.obs[anno_pointdef] == patient_ids[i]
        for j in np.arange(n_subsets):
            ind1 = adata.obs[anno_subsets] == subset_ids[j]
            subset_frac.iloc[i,j] = sum(ind1&ind2)
    subset_frac = subset_frac.apply(lambda x: x/sum(x),axis=1)

    #Fetch patient labels: 
    patient_phenos = pd.DataFrame(index = subset_frac.index, columns=[anno_groupby])
    for i in range(0,len(labels)):
        p_list = adata.obs[anno_pointdef][adata.obs[anno_groupby] == labels[i]].unique().tolist()
        patient_phenos[anno_groupby][p_list] = labels[i]
    patient_phenos['time'] = patient_phenos['time'].astype('category')
    patient_phenos['time'] = patient_phenos['time'].cat.reorder_categories(labels)
    ord_ind = patient_phenos.sort_values('time', ascending=False).index

    subset_frac = subset_frac.loc[ord_ind]

    fig, axes = plt.subplots(1,1, figsize=(9,12))
	
	#Insert the correct color codes: 
    subset_frac.plot.barh(stacked=True, grid=False, legend=False, ax=axes, color= ['#ff0000', '#c0c999', '#9e1a1a', '#721313', '#5c7aff', '#dec1ff', '#8c29ff',
                                 '#fe6776', '#63264a', '#ff99cc', '#ff0080', '#cc33ff', '#a799b7', '#00b3b3', '#004d4d', 
                                 '#fd96A9'])

    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    axes.legend(loc='lower right', bbox_to_anchor=(1.5, 0.5), frameon=True)
    sns.despine(left=True, bottom=True, right=True)
    plt.tight_layout()
    plt.savefig(save_file)



import matplotlib.pyplot as plt
import seaborn as sns

#Create bar-plots per donor-id a/c the function above: make_bar_plots() 
make_bar_plots(adata_hvg, anno_pointdef = 'donor_id_reordered')



#Since vFB2 & vEB2 are coming from single donors respectively, it's better to filter them out from further analysis. 
adata_filter= adata_hvg[~adata_hvg.obs['celltype_v1'].isin(['vFB2', 'vEB2'])]
adata_filter.obs['celltype_v1'].cat.categories
patient_ids = adata_filter.obs['donor_id_reordered'].unique()
print(patient_ids) 


#Copy the make_bar_plots() function: 
def make_bar_plots(adata, anno_groupby = 'time',
                         anno_subsets = 'celltype_v1',
                         anno_pointdef = 'donor_id_reordered',
                         save_file = 'SP14_SP82_SP136_placenta_atlas_comp_filtered_220322.pdf'):
    
    #Get number of categories
    labels = adata.obs[anno_groupby].cat.categories.tolist()
    n_labels = len(labels)
    subset_ids = adata.obs[anno_subsets].cat.categories.tolist()
    n_subsets = len(subset_ids)
    patient_ids = adata.obs[anno_pointdef].cat.categories.tolist()
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

    #Fetch patient labels
    patient_phenos = pd.DataFrame(index = subset_frac.index, columns=[anno_groupby])
    for i in range(0,len(labels)):
        p_list = adata.obs[anno_pointdef][adata.obs[anno_groupby] == labels[i]].unique().tolist()
        patient_phenos[anno_groupby][p_list] = labels[i]
    patient_phenos['time'] = patient_phenos['time'].astype('category')
    patient_phenos['time'] = patient_phenos['time'].cat.reorder_categories(labels)
    ord_ind = patient_phenos.sort_values('time', ascending=False).index

    subset_frac = subset_frac.loc[ord_ind]

    fig, axes = plt.subplots(1,1, figsize=(9,12))

    subset_frac.plot.barh(stacked=True, grid=False, legend=False, ax=axes, color= ['#ff0000', '#c0c999', '#9e1a1a', '#5c7aff', '#dec1ff', '#fe6776',
       '#63264a', '#ff99cc', '#ff0080', '#cc33ff', '#a799b7', '#00b3b3','#004d4d', '#fd96A9'])

    box = axes.get_position()
    axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    axes.legend(loc='lower right', bbox_to_anchor=(1.5, 0.5), frameon=True)
    sns.despine(left=True, bottom=True, right=True)
    plt.tight_layout()
    plt.savefig(save_file)


#Filter the FB2/EB2 out since they're donor-specific clusters & remake the bar-plots (important!):
#Annotations are almost final but certain changes (like APAhi: vCTBpf) were made later. 
make_bar_plots(adata_filter, anno_pointdef = 'donor_id_reordered')


#Split by library (10X V3 vs 10X V2) post filtering vFB2 & vEB2: 
for i in adata_filter.obs['library'].cat.categories:
    print(i) 
    fig= sc.pl.umap(adata_filter[adata_filter.obs['library'] == i], color = 'celltype_v1', return_fig=True, title= i)
    #pdf.savefig(fig)


#Split by donor_id: a qualitative visualization to assess concordance of the batches. 
for i in adata_filter.obs['donor_id_reordered'].cat.categories:
    print(i) 
    fig= sc.pl.umap(adata_filter[adata_filter.obs['donor_id_reordered'] == i], color = 'celltype_v1', return_fig=True, title= i)
    


adata_scvi= adata_hvg.copy()
del adata_scvi.obsm['_scvi_extra_categoricals']
del adata_scvi.obsm['_scvi_extra_continuous']
del adata_scvi.obsm['X_normalized_scVI']

#Write the results:
#Foundation of the final analysis performed on whole placenta dataset: 
results_scvi= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_scVI_harmonization/scvi_integration_1803/SP014_SP082_SP136_villi_scvi_lv15.h5ad"
adata_scvi.write(results_scvi)

##Important note (after September'22): Old path /data/analysis/preeclampsia_2019/.. is changed to /dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/villi_anndata for storage in Eils-HPC cluster. 

#Print the cell type composition of control samples (finalized upon trajectory analysis)
adata_scvi[adata_scvi.obs['disease']== 'C'].obs['celltype_v1'].value_counts()

#Compute cell markers- Naive Bayes method implemented in scVI ("change_mode"): 
#Reference: Deep Generative Models for Detecting Differential Expression in Single Cells (https://www.biorxiv.org/content/biorxiv/early/2019/10/04/794289.full.pdf)  
de_df = vae_lv15.differential_expression(groupby= "celltype_v1")
de_df.to_csv('./scvi_bayes_markers_2103/Placenta_Bayes_corr_celltypev1.csv')
de_df.head()

#Filter vFB2 & vEB2 from the markers data-frame:
options= ['vFB2', 'vEB2']
df_filter = de_df[~de_df['group1'].isin(options)] 
df_filter.head()

#Store the cluster specific features after removing vFB2 & vEB2: 
markers = {}
cats = adata_hvg.obs.celltype_v1.cat.categories

#Store the filtered DEG: from Naive Bayes (used internally for cluster validations)
for i, c in enumerate(cats):
    cid = "{} vs Rest".format(c)
    cell_type_df = df_filter.loc[df_filter.comparison == cid]
	cell_type_df = cell_type_df[cell_type_df.lfc_mean > 1]
	cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3] #cut-off recommended by scVI. 
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.25] #pct= 25% (standard in Seurat)

    markers[c] = cell_type_df
    
pd.concat(markers).to_csv('./scvi_bayes_markers_2103/Placenta_celltypev1_Bayes3_lfc1_pct25_210322.csv')


#Store the filtered DEG upon increasing the "non_zeros_proportion1" (or, percentage of cells expressing a DEG) to 30%. 
markers = {}

for i, c in enumerate(cats):
    cid = "{} vs Rest".format(c)
    cell_type_df = df_filter.loc[df_filter.comparison == cid]
	cell_type_df = cell_type_df[cell_type_df.lfc_mean > 1]
	cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3]
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.30]

    markers[c] = cell_type_df
    
pd.concat(markers).to_csv('./scvi_bayes_markers_2103/Placenta_celltypev1_Bayes3_lfc1_pct30_210322.csv')


#Use a combination of established & computational cell types/states markers: 
cell_markers= ['GREM2', 'ESRRG', 'ERVFRD-1', 'GATA3', 'TACC2', 'PEG10', 'PARP1', 'FRAS1', 'TP63', 'YAP1', 'PBX1', 
               'MKI67', 'TOP2A', 'DIAPH3', 'POLQ','ATAD2', 'CENPE', 'CGA', 'CSH1', 'CSH2', 'TFPI2', 'AFF1', 
               'CYP19A1', 'ADAM12', 'KISS1', 'GH2', 'HLA-G', 'MYCNUT', 'NOTUM', 'DIO2', 'ASCL2', 'LAIR2', 
              'COL27A1', 'MS4A6A', 'CD14', 'MRC1', 'SPP1', 'VSIG4', 'F13A1', 'LYVE1', 'CD36', 'RBPJ', 'MAF', 'CD163', 
               'SIGLEC1', 'HLA-DRB1', 'HLA-DPA1', 'LYZ', 'DPYD', 'CLEC4E', 'CIITA', 'CD48', 'CD44', 'CD83', 'APOE', 'CD74', 
               'CD2', 'CD3G', 'SKAP1', 'THEMIS', 'CD96', 'RUNX3', 'LCK', 'ITK', 'HLA-A', 'HLA-B', 
               'ACTA2', 'GUCY1A2', 'AGTR1', 'AFF2', 'PLA2G5', 'COL5A3', 'HEYL', 'EGFLAM', 'CCDC102B', 
               'EGFL6', 'SOX5', 'COL1A1','COL1A2', 'COL3A1', 'CDH11', 'DCN', 
               'CHSY3', 'LSAMP','POSTN', 'DKK2', 'ROBO2', 'CNTN5', 
               'CD34', 'MEOX2', 'PECAM1', 'LDB2', 'DACH1',
               'HBZ', 'HBA1', 'HBA2', 'EPB41', 'ACSM3', 'TRAK2', 'NARF', 'MICAL2']

sc.pl.dotplot(adata_filter, cell_markers, groupby='celltype_v1', dendrogram=False, 
              color_map="Blues", use_raw=True, standard_scale="var") #Expression scaled from 0-1 for visualization. 

#Convert to full-length dataset: 
ldata_norm_int= adata_filter.raw.to_adata()
sc.pp.normalize_total(ldata_norm_int, target_sum=1e4)
sc.pp.log1p(ldata_norm_int) #Log-transform. 

#Reorder the renamed/existing categories:
adata_filter.obs['celltype_v1_reordered']= adata_filter.obs['celltype_v1'].cat.reorder_categories(['vVCT_prol', 'vVCT', 'APAhi_tropho',
    'vSCTjuv', 'vSCT_1', 'vSCT_2', 'vEVT', 'vFB1', 'vMC', 'vVEC', 'vHBC', 'PAMM', 'vTcell', 'vEB1'])

adata_filter.obs['celltype_v1_reordered'].cat.categories

cell_markers= ['GREM2', 'ESRRG', 'ERVFRD-1', 'ERVV-1', 'PEG10', 'PARP1', 'FRAS1', 'TP63', 'YAP1', 'PBX1', 
               'MKI67', 'TOP2A', 'POLQ','ATAD2', 'CENPE', 'CENPK', 'CGA', 'CSH1', 'CSH2', 'TFPI2', 'ALPP', 
               'CYP19A1', 'ADAM12', 'KISS1', 'GH2', 'HLA-G', 'MYCNUT', 'NOTUM', 'DIO2', 'ASCL2', 'LAIR2', 
              'COL27A1', 'MS4A6A', 'CD14', 'MRC1', 'SPP1', 'VSIG4', 'F13A1', 'LYVE1', 'CD36', 'RBPJ', 'MAF', 'CD163', 
               'SIGLEC1', 'HLA-DRB1', 'HLA-DPA1', 'LYZ', 'DPYD', 'CIITA', 'CD48', 'CD44', 'CD83', 'APOE', 'CD74', 
               'CD2', 'CD3G', 'SKAP1', 'THEMIS', 'CD96', 'RUNX3', 'LCK', 'ITK', 'HLA-A', 'HLA-B', 
               'ACTA2', 'GUCY1A2', 'AGTR1', 'AFF2', 'PLA2G5', 'COL5A3', 'HEYL', 'EGFLAM', 'CCDC102B', 
               'EGFL6', 'SOX5', 'COL1A1','COL1A2', 'COL3A1', 'CDH11', 'DCN',
               'CD34', 'MEOX2', 'PECAM1', 'LDB2', 'DACH1',
               'HBZ', 'HBA1', 'HBA2', 'EPB41', 'ACSM3', 'TRAK2', 'NARF', 'MICAL2']


#Re-plot the cell markers: 
sc.pl.dotplot(adata_filter, cell_markers, groupby='celltype_v1_reordered', dendrogram=False, 
              color_map="Blues", use_raw=True, standard_scale="var", save= '_placenta_markers_2103_v1.pdf')


#Check robustness by fine-tuning the n_neighbors parameter: 
ldata_umap= adata_hvg.copy() 
#Use k=20: 
sc.pp.neighbors(ldata_umap, n_neighbors=20, use_rep= 'X_scVI')
sc.tl.umap(ldata_umap, random_state=0) 
sc.pl.umap(ldata_umap, color= ['celltype_v1'])



#Use k=25. 
ldata_umap01= adata_hvg.copy() 
sc.pp.neighbors(ldata_umap01, n_neighbors=25, use_rep= 'X_scVI')
sc.tl.umap(ldata_umap01, random_state=0) 
sc.pl.umap(ldata_umap01, color= ['celltype_v1'])


adata_scvi= ldata_umap.copy() #Create a copy. 
adata_scvi.raw.X #double check the dimension 
sc.pl.umap(adata_scvi, color= ['celltype_v1'])


#Re-write the results:
results_scvi= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_scVI_harmonization/scvi_integration_1803/SP014_SP082_SP136_villi_scvi_lv15.h5ad"
adata_scvi.write(results_scvi)
del adata_scvi.obsm['X_normalized_scVI']


#Repeat the filtering steps to remove vFB2/vEB2 as explained above: 
adata_filter= ldata_umap[~ldata_umap.obs['celltype_v1'].isin(['vFB2', 'vEB2'])]
sc.pl.umap(adata_filter, color= ['celltype_v1'])


#Reorder the cell type categories. 
#Note that APAhi_tropho was renamed to vCTBpf; vSCT= vSTB, vVCT= vCTB & vEVT= vCCT (proximal cell columns) later. 
adata_filter.obs['celltype_v1_reordered']= adata_filter.obs['celltype_v1'].cat.reorder_categories(['APAhi_tropho', 'vVCT', 'vVCT_prol',
    'vSCTjuv', 'vSCT_1', 'vSCT_2', 'vEVT', 'vFB1', 'vMC', 'vVEC', 'vHBC', 'PAMM', 'vTcell', 'vEB1'])



#Plot senescence markers (after discussing with Dr. Olivia Nonn):
on_markers= ['CDKN1A', 'CDKN2A', 'TP53', 'LAMB1', 'RB1', 'TP53BP1'] #'H2AX'

#Negligible expression of all of them. 
sc.pl.dotplot(adata_filter, on_markers, groupby= ['celltype_v1_reordered'], dendrogram=False, 
              color_map="Blues", use_raw=True, standard_scale="var") 



#check by reading the saved anndata again:
results_scvi= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_scVI_harmonization/scvi_integration_1803/SP014_SP082_SP136_villi_scvi_lv15.h5ad"
adata_villi= sc.read_h5ad(results_scvi)
ldata_norm_int= adata_villi.raw.to_adata() #Convert to full-length matrix. 
sc.pp.normalize_total(ldata_norm_int, target_sum=1e4) #Normalize. 
sc.pp.log1p(ldata_norm_int) #take the logarithm (log1p)


#Subset vSTB subgroups: 
subset= ['vSCT_1', 'vSCT_2', 'vSCTjuv']
adata_sct= ldata_norm_int[ldata_norm_int.obs['celltype_v1'].isin(subset)]
adata_sct.obs['celltype_v1'].value_counts() #print the initial #STB after analysis so far (59214 vSCT1; 13293 vSCT2 & 3392 vSCTjuv)


#Inspect well known vSCT markers: 
sct_gene_list= ['ADAM12', 'SVEP1', 'LINC01949', 'CHODL', 'PAPPA', 'PAPPA2', 'PSG11', 'PSG8', 'PSG1', 
                'KISS1', 'CGA', 'PGF', 'CYP19A1', 'KYNU', 'EBI3', 
                'CYP11A1', 'NECTIN3', 'BACE2', 'COBLL1', 'ATG9B', 'PSG9', 'EZR', 'EGFR', 
                'PSG5', 'DLG5', 'ZNF292', 'GULP1', 'ADGRL3', 'GATA2', 'JUND', 'AC090826.1', 
                'CSH1', 'CSH2', 'TFPI2', 'VIM', 'SPARC', 'DLK1', 'TMSB10', 'ACTB', 'IGFBP3', 'ALPP', 
               'FBLN1', 'MMP11', 'FBN2', 'KRT18', 'KRT19', 'PSG3'] 

sc.set_figure_params(dpi=100)


sc.pl.dotplot(ldata_norm_int[ldata_norm_int.obs['celltype_v1'].isin(subset)], 
              sct_gene_list, groupby='celltype_v1', dendrogram=False, 
              color_map="Blues")



plt.rcParams["figure.figsize"] = (4,4) #re-set figure size.

#Feature plots of vSTB markers; ECM organization genes such as SPARC & DLK1 plotted to check for vSCTjuv. 
sc.pl.umap(ldata_norm_int, color= ['celltype_v1', 'CGA', 'CSH1', 'CSH2', 'CYP19A1', 'TFPI2', 'PSG1',
                                  'PSG2', 'PSG4', 'PSG8', 'SPARC', 'DLK1'])



#Plot the cell annotations (predicted) & Leiden clusters below each other 
sc.pl.umap(ldata_norm_int, color= ['celltype_v1', 'leiden'], legend_loc='on data',
           frameon=False, legend_fontsize=8, legend_fontoutline=2, ncols=1)




#Subset the trophoblast cell groups only: 'vSCT_1', 'vSCT_2', 'vSCTjuv', 'vVCT', 'vVCT_prol', 'vEVT', 'APAhi_tropho'
adata_tropho= adata_villi[adata_villi.obs['celltype_v1'].isin(['vSCT_1', 'vSCT_2', 'vSCTjuv', 'vVCT', 'vVCT_prol', 'vEVT', 'APAhi_tropho'])]

#Further subset the controls only (annotations of which were finalized for trajectory analysis)															  
adata_tropho= adata_tropho[adata_tropho.obs['disease']== 'C']
adata_tropho.obs['leiden'].to_csv('Placenta_scvi_leiden_2403.csv') #Save the Leiden labels of trophoblast. 


#Crosstab to check how the Leidens conform to LR predicted annotations for trophoblast (stored under "vVCT_subclust")
pd.crosstab(adata_tropho.obs['leiden'], adata_tropho.obs['vVCT_subclust']) #L-36 should be VCT or VCTprol 


adata_tropho[adata_tropho.obs.index.duplicated()] #No duplicated index! 



#Filter the Leiden-16.1 out:
adata_filter= adata_villi[adata_villi.obs['leiden_16']!= 'Leiden_16,1']

sc.pl.umap(adata_filter, color= ['celltype_v1'], legend_loc='on data', frameon=False, legend_fontsize=8, legend_fontoutline=2)


#Extract Leiden_16.1 out (steps below are done to refine the vSTBjuv & vCTBp clusters and prevent minor misclassifications)
#A notebook under the same name is added for convenience to understand the steps in details: 
#leiden_16: represents L16 subclustering performed previously on control dataset. 
adata_leiden16= adata_tropho[adata_tropho.obs['leiden_16']== 'Leiden_16,1']
adata_leiden16.n_obs  #547 cells; possibly ambiguous cluster having mixed expression of vSTBjuv, vCTBprol & some vCCT features. 


#Refinement of vSTBjuv cluster: 
juv_leiden16= pd.Series(adata_leiden16.obs['leiden_16']) #Grab the barcodes of leiden 16.1 subcluster. 
adata_villi.obs["leiden_subclusters"]= adata_villi.obs["celltype_v1"] #create a copy. 
adata_villi.obs["leiden_subclusters"] = adata_villi.obs["leiden_subclusters"].astype(str)

#Replace the loci of cells to 36 & 38 leiden: 
adata_villi.obs["leiden_subclusters"].loc[juv_leiden16.index] = juv_leiden16 
adata_villi.obs["leiden_subclusters"].value_counts() 
adata_villi.obs['leiden_subclusters']= adata_villi.obs['leiden_subclusters'].astype('category')

#Here, the Leiden-16.1 (anomalous cluster) can be spotted separately above vSTBjuv & in close contact with vCTB. 
sc.pl.umap(adata_villi, color= ['leiden_subclusters'], legend_loc='on data',
           frameon=False, legend_fontsize=8, legend_fontoutline=2)



#Leiden-16.1 is removed from analysis owing to mixed marker genes: 
adata_filter= adata_villi[adata_villi.obs['leiden_subclusters']!= 'Leiden_16,1']
adata_filter.uns['leiden_subclusters_colors']= adata_filter.uns['celltype_v1_colors']

#UMAP after filtering L16.1 (perfect!)
sc.pl.umap(adata_filter, color= ['leiden_subclusters', 'leiden'], legend_loc='on data', frameon=False, legend_fontsize=8, legend_fontoutline=2)



#L39 has a strong chance of being vHBC_prol (given robust expression of TOP2A/MKI67)
#Interestingly, they are mainly found in first-trimester villi (435 early; 4 pre-term)
pd.crosstab(adata_filter[adata_filter.obs['leiden']== '39'].obs['leiden'], adata_filter.obs['time'])
adata_filter.obs.leiden_subclusters_refined= adata_filter.obs.leiden_subclusters


#Only the early cells of leiden-39 should be mapped to 'vHBCp'
#adata_filter.obs.loc[(adata_filter.obs.time== 'early') & (adata_filter.obs.leiden == '39') & (adata_filter.obs.leiden_subclusters == 'vHBC'), 'leiden_subclusters_refined'] = 'vHBC_prol'
adata_hbc= adata_filter.copy()
adata_hbc= adata_hbc[adata_hbc.obs.leiden == '39'] 
adata_hbc.n_obs #Extracted 439 vHBCp cells. 
 
#Rename the Leiden in HBC subset: 
adata_hbc.obs.leiden= (adata_hbc.obs.leiden.map(lambda x: {'39': 'vHBC_prol'}.get(x, x)).astype("category"))
adata_hbc.obs.leiden= adata_hbc.obs.leiden.astype('str')
hbc_prol= pd.Series(adata_hbc.obs.leiden)

#Add a refined variant of "leiden_subclusters" to store the additional vHBCp. 
adata_filter.obs['leiden_subclusters_refined']= adata_filter.obs['leiden_subclusters'] 
adata_filter.obs['leiden_subclusters_refined']= adata_filter.obs['leiden_subclusters_refined'].astype('string')
#Replace the existing leiden-39 cell barcodes as vHBC_prol:
adata_filter.obs['leiden_subclusters_refined'].loc[hbc_prol.index] = hbc_prol
adata_filter.obs['leiden_subclusters_refined'].value_counts() #Wow! 


#Set new color code for vHBCp: 
adata_filter.uns['leiden_subclusters_refined_colors']= ['#ff0000', '#c0c999', '#9e1a1a', '#721313', '#5c7aff', '#dec1ff',
       '#8c29ff', '#fe6776', '#b32d00', '#63264a', '#ff99cc', '#ff0080', '#cc33ff',
       '#a799b7', '#00b3b3', '#004d4d', '#fd96A9']


#Plot refined annotations after filtering out FB2/EB2: 
sc.pl.umap(adata_filter[~adata_filter.obs['leiden_subclusters_refined'].isin(['vEB2', 'vFB2'])], color= ['leiden_subclusters_refined'], frameon=False, legend_fontsize=8)

#Generate figures & save as PDF:
sc.pl.umap(adata_filter[~adata_filter.obs['leiden_subclusters_refined'].isin(['vEB2', 'vFB2'])], 
           color= ['leiden_subclusters_refined'], 
           frameon=False, legend_fontsize=8, save= '_placenta_umap_final_240322.pdf')

adata_filter02= adata_filter[~adata_filter.obs['leiden_subclusters_refined'].isin(['vEB2', 'vFB2'])]

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
#Split UMAP by 10X V3 & V2: 
for i in adata_filter02.obs['library'].cat.categories:
    print(i) 
    fig= sc.pl.umap(adata_filter02[adata_filter02.obs['library'] == i], color = 'leiden_subclusters_refined', return_fig=True, title= i)
    #pdf.savefig(fig)


#Map the whole L28 as 'vMC': 
#The vSCT2 cells within L28 is a result of misclassifion. 
adata_filter02.obs.loc[(adata_filter02.obs.leiden == '28'), 'leiden_subclusters_refined'] = 'vMC'
adata_filter02.obs['leiden_subclusters_refined'].value_counts() 
adata_filter02[adata_filter02.obs.leiden == '3'].obs['leiden_subclusters_refined'].value_counts()


#Remove the APAhi_tropho (vCTBpf) in VCTprol (vCTBp) cluster (aka leiden-3)
adata_filter02.obs.leiden_subclusters_refined02= adata_filter02.obs.leiden_subclusters_refined
adata_filter02.obs.leiden_subclusters_refined02= adata_filter02.obs.leiden_subclusters_refined02.astype('object')
adata_filter02.obs.loc[(adata_filter02.obs.leiden == '3') & (adata_filter02.obs.leiden_subclusters_refined == 'APAhi_tropho'),
                       'leiden_subclusters_refined02'] = 'vVCT_prol'

adata_filter02[adata_filter02.obs.leiden == '3'].obs['leiden_subclusters_refined'].value_counts() 
adata_filter02[adata_filter02.obs.leiden == '3'].obs['leiden_subclusters_refined02'].value_counts() #worked!  


sc.settings.set_figure_params(dpi=120)
plt.rcParams["figure.figsize"] = (3,3)
#Generate figures for the manuscript:
sc.pl.umap(adata_filter02, color= ['leiden_subclusters_refined'],  frameon=False, legend_fontsize=8)

#The workaround & annotation fixes after line-535 is coupled in a notebook & shared under: Placenta_analysis_manuscript_2022.ipynb 




