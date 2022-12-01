#!/usr/bin/env python
#Analysis performed by: Olivia Debnath
#Computational Oncology/BIH (Dr. Naveed Ishaque & Prof. Dr. Roland Eils)
import numpy as np
import pandas as pd
import scanpy as sc

#Import scvi; the global seed is automatically set to 0 for reproducibility. 
import scvi

#Read Cellbender filtered object- refer to the QC notebook to see how the H5AD anndata was built. 
#Important note (after September'22): Old path /data/analysis/preeclampsia_2019/.. is changed to /dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/villi_anndata for storage in Eils-HPC cluster. 
latest= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_anndata/Placenta_cellbender_normalized_matrix_270122.h5ad"

#Read the anndata:
ldata= sc.read_h5ad(latest)


#Scan for NaN or infinity values. 
X = ldata.X
X = np.log1p(X)
X = np.expm1(X)
mean = X.mean(axis=0).A1
mean[mean == 0] = 1e-12
mean = np.log1p(mean)
np.any(mean == np.inf) #No Inf/NaN exists in our matrix. 


#Compute top 4000 HVG using donor_ID as a batch key to prevent high sample specific gene variability. 
#Features variable in at least two donors were considered for downstream analysis. 
#Note, the steps below is for the cell typing for control villi & don't represent the final analysis/UMAP in the manuscript. 
sc.pp.highly_variable_genes(ldata, n_top_genes= 4000, batch_key= "donor_id", subset=True)


#Critical step after HVG computation: set up the anndata to run the VAE on. 
#We will start with a simple model training (Case-I) without encoding the extra covariates.   
scvi.data.setup_anndata(ldata, layer="counts", batch_key= "donor_id")

#Case-I: Model training using 15 latent variables. 
#n_latent: Dimensionality of the latent space
#dropout_rate= 0.1 (default; Dropout rate for neural networks)
#gene_likelihood: zinb;  Zero-inflated negative binomial distribution. 
#n_layers: Number of hidden layers used for encoder and decoder NNs.
#VAE parameters (default): https://docs.scvi-tools.org/en/0.9.1/api/reference/scvi.module.VAE.html#scvi.module.VAE 

#n_latent= 10, 15, 20, 30, 40 showed similar data structures as visualized by UMAP. 
#The following model isn't finally used. Check below for the data harmonization after encoding different covariates. 
vae_lv15 = scvi.model.SCVI(ldata, n_layers=2, n_latent=15)
vae_lv15.train()

#Print model attributes 
vae_lv15 


#Get the latent representation from model training to learn the nearest neighbor graph. 
ldata.obsm["X_scVI"] = vae_lv15.get_latent_representation()
ldata.obsm["X_normalized_scVI"] = vae_lv15.get_normalized_expression()
#Finally, we can compute the KNN graph and visualize it with UMAP:
#Run with default parameters but fix a random seed. 
sc.pp.neighbors(ldata, use_rep="X_scVI")
sc.tl.leiden(ldata, resolution= 1.5, random_state= 0)
sc.tl.umap(ldata, random_state= 0)


#Read the already (basic) predictions on control data. Performed using Celltypist on our old dataset (SP014; from 2021 submission). 
#Reference notebook: Placenta_celltypist_control_predictions_280122_S1.ipynb  
results_annot= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/placenta_celltypist_model/Placenta_atlas_celltypist_annotated_2701.h5ad" 
adata_annot= sc.read_h5ad(results_annot)

#Transfer the celltypist predicted_labels & majority voting: 
ldata.obs['celltypist_majority_voting']= adata_annot.obs['majority_voting']
ldata.obs['celltypist_predicted_labels']= adata_annot.obs['predicted_labels']



#Set the color codes per cluster: 
ldata.uns['celltypist_majority_voting_colors']= ["#048bA8", "#c0c999", "#9e1a1a", "#5c7aff", "#dec1ff", 
                                                 "#fe6776", "#63264a", "#ff99cc", "#ff0080", "#ff0000", 
                                                 "#a799b7", "#00b3b3", "#004d4d", "#fd96A9", "#cc33ff"]
#Switch to refined colors: 
ldata.uns['celltypist_predicted_labels_colors']= ["#048bA8", "#c0c999", "#9e1a1a", "#5c7aff", "#dec1ff", 
                                                 "#fe6776", "#63264a", "#ff99cc", "#ff0080", "#ff0000", 
                                                 "#a799b7", "#00b3b3", "#004d4d", "#fd96A9", "#cc33ff"]




#Set figure parameters: 
sc.settings.set_figure_params(dpi=80)
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (4,4)
sc.pl.umap(ldata, color=["donor_id", "celltypist_predicted_labels", "leiden"], ncols=1)




import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

#Split the UMAPs: old (SP014) vs new samples (SP082) to investigate batch-effects. 
pdf = matplotlib.backends.backend_pdf.PdfPages("UMAP_batchID_SP014_SP082_LV15_29012022.pdf")

for i in ldata.obs['time'].cat.categories:
    print(i) 
    fig= sc.pl.umap(ldata[ldata.obs['time'] == i], color = 'celltypist_predicted_labels', return_fig=True, title= i+' samples')
    pdf.savefig(fig)
    
pdf.close()



#Case-II: analysis used for control villi/placenta annotations: 
#Rerun the scvi integration: encode categorical & continuous covariates.
#For few early samples (SP014), BMI/maternal age information were unavailable during January-March'2022.  
scvi.data.setup_anndata(ldata, layer="counts", batch_key= "donor_id",
                       categorical_covariate_keys= ["library", "time", "gestational_weeks"],
                       continuous_covariate_keys=["n_genes_by_counts", "total_counts", "pct_counts_MT_genes", "XIST-counts"])



#Rerun the VAE with LV=15 & 2-layers architecture after encoding the extra covariates: 
vae_lv15 = scvi.model.SCVI(ldata, n_layers=2, n_latent=15)
vae_lv15.train()



#Latent representation: after correcting for batch specific effects (batch_key) plus 3 categorical as well as 4 continuous covariates. 
ldata.obsm["X_scVI"] = vae_lv15.get_latent_representation()
ldata.obsm["X_normalized_scVI"] = vae_lv15.get_normalized_expression()


#Finally, we learn the KNN graph using the latent representation and visualize it using UMAP:
sc.pp.neighbors(ldata, use_rep="X_scVI")
sc.tl.leiden(ldata, resolution= 1.5)
sc.tl.umap(ldata, random_state=0) #Next, increase the maxiter=200, 500 & 1000 to ensure convergence of UMAP algorithm. 

#Visualize by procurement (Graz or Oslo); library (10X V2 vs 10X V3) & initial celltypist predicted labels. 
sc.pl.umap(ldata, color=["procurement", "library", "celltypist_predicted_labels"], ncols=1)


#Split the UMAP by time/disease: 
for i in ldata.obs['time'].cat.categories:
    print(i) 
    fig= sc.pl.umap(ldata[ldata.obs['time'] == i], color = 'celltypist_predicted_labels', return_fig=True, title= i+' samples')
    #pdf.savefig(fig)
    


#Get the reconstruction error for the VAE (after training). 
vae_lv15.get_reconstruction_error()

#Plotting the likelihood change across the 500 epochs of training: blue for training error and orange for testing error.
elbo_train_set = vae_lv15.trainer.history["elbo_train_set"]
elbo_test_set = vae_lv15.trainer.history["elbo_test_set"]
x = np.linspace(0, 400, (len(elbo_train_set)))
plt.plot(x, elbo_train_set, label="train")
plt.plot(x, elbo_test_set, label="test")
plt.ylim(1500, 3000)
plt.legend()



#Visualize the Leiden mappings & compare (qualitatively) with Celltypist predicted labels: 
#Note that all the Leidens were manually inspected using established cell markers (literature plus computational analysis)
#before finalizing annotations. 
sc.pl.umap(ldata, color=["celltypist_predicted_labels", "leiden"], legend_loc= "on data", 
          legend_fontweight= 'normal', legend_fontsize= 'x-small')




#Remove anomalous/spurious Leiden clusters mapping solely to donors/library/technical doublets based:
ldata_umap= ldata[~ldata.obs['leiden'].isin(['29', '28', '23', '32', '31', '21', '30']),:]

#Recompute UMAP after removal of anomalous Leidens: 
sc.tl.umap(ldata_umap, random_state=0) #Next, increase the maxiter=1000. 

#Re-visualize & inspect the Leidens: 
sc.pl.umap(ldata_umap, color=["celltypist_predicted_labels", "leiden"], legend_loc= "on data", 
          legend_fontweight= 'normal', legend_fontsize= 'x-small')


#Compute Differentially Expressed features using scVI's Naive Bayes factor:
#Advantage: query the trained VAE (batch/covariates corrected) for finding cell state specific featutes. 
de_df = vae_lv15.differential_expression(groupby= "celltypist_predicted_labels")
de_df.to_csv('Placenta_atlas_scVI_vae_DE_Bayes_280122.csv') 


#We now extract top markers for each cluster using the DE results to investigate clusters.
markers = {}
cats = ldata_umap.obs.celltypist_predicted_labels.cat.categories
for i, c in enumerate(cats):
    cid = "{} vs Rest".format(c)
    cell_type_df = de_df.loc[de_df.comparison == cid]

    cell_type_df = cell_type_df[cell_type_df.lfc_mean > 0] #lfc_mean can be increased to scan for specific features. 

    cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3] #Taken the recommended bayes_factor cut-off. 
    cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1] #The non_zeros_proportion1 or pct can be increased to 25-30% to scan for robust features. 

    markers[c] = cell_type_df.index.tolist()[:6] #Store the top-6 markers per cluster to visualize via dot-plot below. 


sc.tl.dendrogram(ldata_umap, groupby= "celltypist_predicted_labels", use_rep="X_scVI")

#Scale the features across group on a scale of 0-1 for dot-plot visualization: 
sc.pl.dotplot(ldata_umap, markers, groupby='celltypist_predicted_labels', dendrogram=True,
    color_map="Blues", use_raw=True, standard_scale="var")



#Compute DEG per Leiden groups using scVI's Bayes factor:
#Markers were inspected to finalize annotations later. 
de_leiden = vae_lv15.differential_expression(groupby= "leiden")
de_df.to_csv('Placenta_atlas_scVI_vae_DE_Bayes_leiden_280122.csv') 



#Increase the maxiter to 1000 to assure better convergence. 
ldata_umap01= ldata_umap.copy()
sc.tl.umap(ldata_umap01, random_state=0, maxiter= 1000) #also tried, maxiter= 100, 200 & 500.  
sc.pl.umap(ldata_umap01, color=["celltypist_predicted_labels", "leiden"], legend_loc= "on data", 
          legend_fontweight= 'normal', legend_fontsize= 'x-small')


#Clearly, a part of peripheral cluster of EB is contaminated with vSCT2 markers (results in Leiden-11 to merge with Leiden-6)  
#when the maxiter is increased to 1000. 
#On a good note, we see transition paths from VCT to SCT via: SCTjuvenile & APAhi (formerly TSCs). 


#Re-split the UMAP by time/disease: to investigate batch-effects: 
for i in ldata_umap01.obs['time'].cat.categories:
    print(i) 
    fig= sc.pl.umap(ldata_umap01[ldata_umap01.obs['time'] == i], color = 'celltypist_predicted_labels', return_fig=True, title= i+' samples')
    #pdf.savefig(fig)
    


#Look at the source of contaminating EB/SCT2 cells: is it specific sample_id?! 
#New early samples are all 10X_V3. 
for i in ldata_umap01.obs['library'].cat.categories:
    print(i) 
    fig= sc.pl.umap(ldata_umap01[ldata_umap01.obs['library'] == i], color = 'celltypist_predicted_labels', return_fig=True, title= i+' samples')
    #pdf.savefig(fig)
    



ldata_umap01.raw= ldata_umap01  #freeze the state in `.raw` to keep the full data matrix safe. 

#Normalize & log-transform: 
sc.pp.normalize_total(ldata_umap01, target_sum=10e4)
sc.pp.log1p(ldata_umap01)


#Plot the well-known STB genes to check how they express. 
sct_markers= ['PAPPA', 'PAPPA2', 'ADAM12', 'KISS1', 'CYP19A1', 'CGA', 
              'TENT5A', 'ADGRL3', 'GATA2', 'BACE2', 'ATG9B', 'TFPI2', 'PSG1', 'PSG3', 'PSG4', 'PSG6', 'PSG8', 
             'CSH1', 'CSH2', 'EEF1A1', 'DLK1']

sc.pl.umap(ldata_umap01, color= sct_markers) 



#Compute DEG using scVI's Bayes factor: repeat this step. 
de_leiden = vae_lv15.differential_expression(groupby= "leiden")
de_df.to_csv('Placenta_atlas_scVI_vae_DE_Bayes_leiden_310122.csv') #to be repeated. 


#Save the trained VAE model (LV15 based) for future reference: 
vae_lv15.save("saved_model/")
#model = SCVI.load("saved_model/")


#Inspect how the Leidens are conforming to predicted cell annotations (using celltypist LR classifier)
df= pd.crosstab(ldata_umap01.obs['celltypist_predicted_labels'], ldata_umap01.obs['leiden'])

#Compute a confusion matrix
conf_mat = df/df.sum(axis=1).values[:, np.newaxis]

plt.figure(figsize=(8, 8)) #Set the figure size 
_ = plt.pcolor(conf_mat)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("Predicted: LR classifier supervised")
plt.ylabel("Observed: Leiden unsupervised")




#Leiden-25 is an admixture of VCT/VCTp as well as SCT1. But this is not robust for CSH1/CSH2.
#So, most probably it's not SCTjuv. Investigated separated using trajectory & further analysis.
#sc.pl.umap(ldata_norm, color= ['PAPPA', 'PAPPA2', 'ADAM12', 'KISS1', 'CYP19A1', 'CGA', 
              #'KANK1', 'PEG10', 'TEAD1', 'YAP1', 'PBX1', 'TP63']) 



del ldata_umap.uns['_scvi'] #delete the slot to save anndata. 



#Save the anndata for future reference & finalzing the control villi annotations.
#Note that, this is the output anndata from the basic integration & after filtering of spurious Leiden clusters. 
latest= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_anndata/Placenta_cellbender_LV15_basic.h5ad"
ldata_umap.write(latest)

#More refinements on the control is performed here: Placenta_scVI_controls_celltyping_LV15_310122.ipynb/html or, Placenta_scVI_controls_exploratory_celltyping_LV15_S3.ipynb
#Important note: Old path /data/analysis/preeclampsia_2019/.. is changed to /dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/villi_anndata for storage in Eils-HPC cluster. 
