#!/usr/bin/env python
# coding: utf-8
#Author: Olivia Debnath
#Computational Oncology/BIH (Dr. Naveed Ishaque & Prof. Dr. Roland Eils)


import numpy as np
import pandas as pd
import scanpy as sc

#Note that every sample shown here was processed using 10X V2 library chemistry. 
#Read the Cellbender corrected data for SP014 early decidua: 17-025 decidua (early). 
adata_17025v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_17-025-Decidua_cellbender_filtered.h5')
adata_17025v.var_names_make_unique() 
adata_17025v.obs['donor_id']= 'Donor-17025-decidua'
adata_17025v.obs['time']= 'early'
adata_17025v.obs['disease']= 'C'
adata_17025v.obs['tissue']= 'decidua'
adata_17025v.obs['smoking']= 'No'
adata_17025v.obs['cohort']= 'SP014'
adata_17025v.obs['library']= '10X v2'
adata_17025v.obs['procurement']= 'Graz'
adata_17025v.obs['placental_volume']= 'NA'
adata_17025v.obs['gestational_weeks']= 7 
adata_17025v.obs['gestational_days']= 54 
adata_17025v.obs['maternal_BMI']= 'NA' #BMI unknown at this point. 
adata_17025v.obs['maternal_age']= 'NA'  #age unknown at this point.

#Read the Cellbender corrected SEKR2 decidua (early)
adata_sekr2v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_SEKR2-Decidua_cellbender_filtered.h5')
adata_sekr2v.var_names_make_unique() 
adata_sekr2v.obs['donor_id']= 'Donor-SEKR2-decidua'
adata_sekr2v.obs['time']= 'early'
adata_sekr2v.obs['disease']= 'C'
adata_sekr2v.obs['tissue']= 'decidua'
adata_sekr2v.obs['smoking']= 'No'
adata_sekr2v.obs['cohort']= 'SP014'
adata_sekr2v.obs['library']= '10X v2'
adata_sekr2v.obs['procurement']= 'Graz'
adata_sekr2v.obs['placental_volume']= 'NA'
adata_sekr2v.obs['gestational_weeks']= 7 #approx 
adata_sekr2v.obs['gestational_days']= 49 
adata_sekr2v.obs['maternal_BMI']= 'NA' #BMI unknown at this point. 
adata_sekr2v.obs['maternal_age']= 'NA' #age unknown at this point.

#Read the Cellbender corrected SOZE1-decidua (early): 
adata_soze1v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_SOZE1-Decidua_cellbender_filtered.h5')
adata_soze1v.var_names_make_unique() 
adata_soze1v.obs['donor_id']= 'Donor-SOZE1-decidua'
adata_soze1v.obs['time']= 'early'
adata_soze1v.obs['disease']= 'C'
adata_soze1v.obs['tissue']= 'decidua'
adata_soze1v.obs['smoking']= 'No'
adata_soze1v.obs['cohort']= 'SP014'
adata_soze1v.obs['library']= '10X v2'
adata_soze1v.obs['procurement']= 'Graz'
adata_soze1v.obs['placental_volume']= 'NA'
adata_soze1v.obs['gestational_weeks']= 9 #approx 
adata_soze1v.obs['gestational_days']= 63 
adata_soze1v.obs['maternal_BMI']= 'NA' #BMI unknown at this point. 
adata_soze1v.obs['maternal_age']= 'NA' #age unknown at this point. 

#Read SP014 late term control sample: 327-decidua. 
adata_327v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_327-Decidua_cellbender_filtered.h5')
adata_327v.var_names_make_unique() 
adata_327v.obs['donor_id']= 'Donor-327-decidua'
adata_327v.obs['time']= 'late_term'
adata_327v.obs['disease']= 'C'
adata_327v.obs['tissue']= 'decidua'
adata_327v.obs['smoking']= 'No'
adata_327v.obs['cohort']= 'SP014'
adata_327v.obs['library']= '10X v2'
adata_327v.obs['procurement']= 'Oslo'
adata_327v.obs['placental_volume']= 'NA'
adata_327v.obs['gestational_weeks']= 39
adata_327v.obs['gestational_days']= 275
adata_327v.obs['maternal_BMI']= 22.68  
adata_327v.obs['maternal_age']= 35

#Read SP014 late term control sample: 328-decidua. 
adata_328v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_328-Decidua_cellbender_filtered.h5')
adata_328v.var_names_make_unique() 
adata_328v.obs['donor_id']= 'Donor-328-decidua'
adata_328v.obs['time']= 'late_term'
adata_328v.obs['disease']= 'C'
adata_328v.obs['tissue']= 'decidua'
adata_328v.obs['smoking']= 'No'
adata_328v.obs['cohort']= 'SP014'
adata_328v.obs['library']= '10X v2'
adata_328v.obs['procurement']= 'Oslo'
adata_328v.obs['placental_volume']= 'NA'
adata_328v.obs['gestational_weeks']= 40 
adata_328v.obs['gestational_days']= 278
adata_328v.obs['maternal_BMI']= 25.78 
adata_328v.obs['maternal_age']= 38


#Read SP014 late term control sample: 372-decidua. 
adata_372v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_372-Decidua_cellbender_filtered.h5')
adata_372v.var_names_make_unique() 
adata_372v.obs['donor_id']= 'Donor-372-decidua'
adata_372v.obs['time']= 'late_term'
adata_372v.obs['disease']= 'C'
adata_372v.obs['tissue']= 'decidua'
adata_372v.obs['smoking']= 'No'
adata_372v.obs['cohort']= 'SP014'
adata_372v.obs['library']= '10X v2'
adata_372v.obs['procurement']= 'Oslo'
adata_372v.obs['placental_volume']= 'NA'
adata_372v.obs['gestational_weeks']= 40 
adata_372v.obs['gestational_days']= 283
adata_372v.obs['maternal_BMI']= 23.05 
adata_372v.obs['maternal_age']= 37

#Read the SP014 PE data: 389-decidua. 
adata_389v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_389-Decidua_cellbender_filtered.h5')
adata_389v.var_names_make_unique() 
adata_389v.obs['donor_id']= 'Donor-389-decidua'
adata_389v.obs['time']= 'late_preterm'
adata_389v.obs['disease']= 'PE'
adata_389v.obs['tissue']= 'decidua'
adata_389v.obs['smoking']= 'No'
adata_389v.obs['cohort']= 'SP014'
adata_389v.obs['library']= '10X v2'
adata_389v.obs['procurement']= 'Oslo'
adata_389v.obs['placental_volume']= 'NA' #placental wt= 220 
adata_389v.obs['gestational_weeks']= 27 
adata_389v.obs['gestational_days']= 192
adata_389v.obs['maternal_BMI']= 23.67 
adata_389v.obs['maternal_age']= 26

#Read the SP014 PE data: 419-decidua. 
adata_419v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_419-Decidua_cellbender_filtered.h5')
adata_419v.var_names_make_unique() 
adata_419v.obs['donor_id']= 'Donor-419-decidua'
adata_419v.obs['time']= 'late_preterm'
adata_419v.obs['disease']= 'PE'
adata_419v.obs['tissue']= 'decidua'
adata_419v.obs['smoking']= 'No'
adata_419v.obs['cohort']= 'SP014'
adata_419v.obs['library']= '10X v2'
adata_419v.obs['procurement']= 'Oslo'
adata_419v.obs['placental_volume']= 'NA' #placental wt= 280. 
adata_419v.obs['gestational_weeks']= 34 
adata_419v.obs['gestational_days']= 236
adata_419v.obs['maternal_BMI']= 21.88
adata_419v.obs['maternal_age']= 33

#Read the SP014 PE data: 274-decidua.  
adata_274v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_274-Decidua_cellbender_filtered.h5')
adata_274v.var_names_make_unique() 
adata_274v.obs['donor_id']= 'Donor-274-decidua'
adata_274v.obs['time']= 'late_preterm'
adata_274v.obs['disease']= 'PE'
adata_274v.obs['tissue']= 'decidua'
adata_274v.obs['smoking']= 'No'
adata_274v.obs['cohort']= 'SP014'
adata_274v.obs['library']= '10X v2'
adata_274v.obs['procurement']= 'Oslo'
adata_274v.obs['placental_volume']= 'NA' #placental wt= 280. 
adata_274v.obs['gestational_weeks']= 32 
adata_274v.obs['gestational_days']= 227
adata_274v.obs['maternal_BMI']= 23.18 
adata_274v.obs['maternal_age']= 43

#Merge all early, late & PE decidua into one object:
adata= adata_17025v.concatenate(adata_sekr2v, adata_soze1v, adata_327v, adata_328v, adata_372v, adata_389v, adata_419v, adata_274v)
print(pd.crosstab(adata.obs['time'], adata.obs['library'])) 

#mitochondrial genes
adata.var['MT_genes'] = adata.var_names.str.startswith('MT-') 
#ribosomal genes
adata.var['Ribo_genes'] = adata.var_names.str.startswith(("RPS","RPL"))
#hemoglobin genes.
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

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (6,6)
#Before filtering, visualize the QC per sample: 
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'donor_id', rotation= 90)
#Group by disease/time/library points:
plt.rcParams["figure.figsize"] = (4,6)

#Before filtering: 
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'disease', rotation= 90)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'time', rotation= 90)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'library', rotation= 90)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'procurement', rotation= 90)

#Filtering: A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. 
#Here we will only consider cells with at least 200 detected genes and genes need to be expressed in at least 3 cells. 
#Please note that those values are highly dependent on the library preparation method used & isolation process.
adata_new= adata.copy()
#Perform all filtering on adata_new: 
sc.pp.filter_cells(adata_new, min_genes=200)
sc.pp.filter_genes(adata_new, min_cells=3)
print(adata_new.n_obs, adata_new.n_vars)
sc.pl.highest_expr_genes(adata_new, n_top=30) #MALAT1 is an ubiquitous gene. 

#Filter for percent mito: keep 5% cut-off. (stick to this)
#The old data will be reanalyzed with 5% cut-off. 
adata_mito_filtered = adata_new[adata_new.obs['pct_counts_MT_genes'] < 5, :]
#After filtering: 
sc.pl.highest_expr_genes(adata_mito_filtered, n_top=30) #MALAT1 is an ubiquitous gene. 

#After filtering: important QC plot. 
sc.pl.violin(adata_mito_filtered, ['pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'donor_id', rotation= 90)
sc.pl.violin(adata_mito_filtered, ['pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'disease', rotation= 90)
sc.pl.violin(adata_mito_filtered, ['pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'time', rotation= 90)
sc.pl.violin(adata_mito_filtered, ['pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'procurement', rotation= 90)

#Summarize the decidua QC:: groupby disease. 
sc.pl.violin(adata_mito_filtered, ['n_genes_by_counts', 'total_counts', 'pct_counts_MT_genes','pct_counts_Ribo_genes'],
             jitter=0, groupby = 'disease', rotation= 90, save= '_decidua_disease_QC_270122.pdf')

#Color inputs must be from either .obs or .var, so add in XIST expression to obs.
adata_mito_filtered.obs["XIST-counts"] = adata_mito_filtered.X[:,adata_mito_filtered.var_names.str.match('XIST')].toarray()
sc.pl.violin(adata_mito_filtered, ['XIST-counts'],jitter=0, groupby = 'donor_id', rotation= 90)

#Make scatter plots visualizing how total counts relate to #genes per cell. 
#Vienna recommendation: Don't blindly remove cells by "total counts" => might remove the tetraploid EVTs. 
#No gold standard cut-off for the placenta. 
sc.pl.scatter(adata_mito_filtered, x='total_counts', y='pct_counts_MT_genes')
sc.pl.scatter(adata_mito_filtered, x='total_counts', y='n_genes_by_counts')

#Save the filtered raw matrix: 
results= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_anndata/Decidua_cellbender_raw_matrix_270122.h5ad"
adata_mito_filtered.write(results)

#Calculate cell-cycle scores:
#We, here perform cell cycle scoring. To score a gene list, the algorithm calculates the difference of mean expression of the given list and the mean expression of reference genes. To build the reference, the function randomly chooses a bunch of genes matching the distribution of the expression of the given list. Cell cycle scoring adds three slots in data, a score for S phase, a score for G2M phase and the predicted cell cycle phase.
#First read the file with cell cycle genes, from Regev lab and split into S and G2M phase genes. Cell cycle genes were retrieved from the scanpy_usage github site via web browser at RegevLab Github repo
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
#Scale data: please don't do this step now. scVI requires raw data. 
sc.tl.score_genes_cell_cycle(adata_final, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pl.violin(adata_final, ['S_score', 'G2M_score'],jitter=0, groupby = 'donor_id', rotation=90)
#Make a count layer to securely store the raw data to proceed with scVI harmonization: 
adata_final.layers["counts"] = adata_final.raw.X.copy()

#Re-write the concatenated Cellbender processed, QC filtered & log-normalized object for SP014.
#Raw data stored in count slot for scVI.  
#Latest directory on Eils-HPC (after Sept'22 migration): /dh-projects/preeclampsia_2022/analysis/placenta_atlas_2022/decidua_anndata/Decidua_cellbender_normalized_matrix_270122.h5ad 
latest= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/decidua_anndata/Decidua_cellbender_normalized_matrix_270122.h5ad"
adata_final.write(latest) 

#Note that, same parameters for cell/gene filtering is applied for 10X V3 processed SP136 samples (n=3). 
#Processing of SP136 samples & their subsequent integration with SP014 samples (10X V2 processed) post Cellbender correction is shown here: SP136_SP014_cellbender_decidua_celltyping_2022.py 






