#!/usr/bin/env python
# coding: utf-8
#Author: Olivia Debnath 

import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.set_figure_params(dpi=80)

#Metadata as available in January'2022 (verified from scRNAseq_metadata_early_controls.xlsx; obtained from Dr. Florian Herse on 13.11.2021): 
#Read SP082_011: load Cellbender filtered matrices (after removing background). 
#Latest donor_id: 18-082-01-v
adata_d11 = sc.read_10x_h5('./SP082_011/outs/cellbender_filtered.h5')
adata_d11.var_names_make_unique() 
adata_d11.obs['donor_id']= 'Donor11'
adata_d11.obs['time']= 'early'
adata_d11.obs['disease']= 'C'
adata_d11.obs['tissue']= 'villi'
adata_d11.obs['smoking']= 'No'
adata_d11.obs['placental_volume']= 4.2
adata_d11.obs['gestational_weeks']= 10 
adata_d11.obs['gestational_days']= 74
adata_d11.obs['maternal_BMI']= 22.77
adata_d11.obs['maternal_age']= 28

#Read SP082_012 (18-108-02-v): 
adata_d12 = sc.read_10x_h5('./SP082_012/outs/cellbender_filtered.h5')
adata_d12.var_names_make_unique() 
adata_d12.obs['donor_id']= 'Donor12'
adata_d12.obs['time']= 'early'
adata_d12.obs['disease']= 'C'
adata_d12.obs['tissue']= 'villi'
adata_d12.obs['smoking']= 'No'
adata_d12.obs['placental_volume']= 4
adata_d12.obs['gestational_weeks']= 10 
adata_d12.obs['gestational_days']= 76
adata_d12.obs['maternal_BMI']= 21.01
adata_d12.obs['maternal_age']= 25

#Received in Berlin-Buch Aug 2021:  
#Read P1200_SP082_034 (18-032-v): 
adata_d34 = sc.read_10x_h5('./P1200_SP082_034/outs/cellbender_filtered.h5')
adata_d34.var_names_make_unique() 
adata_d34.obs['donor_id']= 'Donor34'
adata_d34.obs['time']= 'early'
adata_d34.obs['disease']= 'C'
adata_d34.obs['tissue']= 'villi'
adata_d34.obs['smoking']= 'No'
adata_d34.obs['placental_volume']= 'NA'
adata_d34.obs['gestational_weeks']= 8  
adata_d34.obs['gestational_days']= 62 
adata_d34.obs['maternal_BMI']= 17.36
adata_d34.obs['maternal_age']= 20

#Read P1200_SP082_039 (18-098-v): 
adata_d39 = sc.read_10x_h5('./P1200_SP082_039/outs/cellbender_filtered.h5')
adata_d39.var_names_make_unique() 
adata_d39.obs['donor_id']= 'Donor39'
adata_d39.obs['time']= 'early'
adata_d39.obs['disease']= 'C'
adata_d39.obs['tissue']= 'villi'
adata_d39.obs['smoking']= 'No'
adata_d39.obs['placental_volume']= 'NA'
adata_d39.obs['gestational_weeks']= 5 
adata_d39.obs['gestational_days']= 35 
adata_d39.obs['maternal_BMI']= 21.01
adata_d39.obs['maternal_age']= 28
adata_d39

#Read P1200_SP082_041 (20-027-v): 
adata_d41 = sc.read_10x_h5('./P1200_SP082_041/outs/cellbender_filtered.h5')
adata_d41.var_names_make_unique() 
adata_d41.obs['donor_id']= 'Donor41'
adata_d41.obs['time']= 'early'
adata_d41.obs['disease']= 'C'
adata_d41.obs['tissue']= 'villi'
adata_d41.obs['smoking']= 'No'
adata_d41.obs['placental_volume']= 'NA'
adata_d41.obs['gestational_weeks']= 6 
adata_d41.obs['gestational_days']= 42 
adata_d41.obs['maternal_BMI']= 19.72
adata_d41.obs['maternal_age']= 32

#Read P1200_SP082_037 (18-017-v): 
adata_d37 = sc.read_10x_h5('./P1200_SP082_037/outs/cellbender_filtered.h5')
adata_d37.var_names_make_unique() 
adata_d37.obs['donor_id']= 'Donor37'
adata_d37.obs['time']= 'early'
adata_d37.obs['disease']= 'C'
adata_d37.obs['tissue']= 'villi'
adata_d37.obs['smoking']= 'No'
adata_d37.obs['placental_volume']= 'NA'
adata_d37.obs['gestational_weeks']= 8 
adata_d37.obs['gestational_days']= 58 
adata_d37.obs['maternal_BMI']= 25.18
adata_d37.obs['maternal_age']= 35

#Read P1200_SP082_038 (18-033-v): 
adata_d38 = sc.read_10x_h5('./P1200_SP082_038/outs/cellbender_filtered.h5')
adata_d38.var_names_make_unique() 
adata_d38.obs['donor_id']= 'Donor38'
adata_d38.obs['time']= 'early'
adata_d38.obs['disease']= 'C'
adata_d38.obs['tissue']= 'villi'
adata_d38.obs['smoking']= 'No'
adata_d38.obs['placental_volume']= 'NA'
adata_d38.obs['gestational_weeks']= 9 
adata_d38.obs['gestational_days']= 64 
adata_d38.obs['maternal_BMI']= 21.01 
adata_d38.obs['maternal_age']= 32

#Read the Cellbender corrected data for SP014 early: 17-025-v. 
adata_17025v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_17-025-Villi_cellbender_filtered.h5')
adata_17025v.var_names_make_unique() 
adata_17025v.obs['donor_id']= 'Donor-17025'
adata_17025v.obs['time']= 'early'
adata_17025v.obs['disease']= 'C'
adata_17025v.obs['tissue']= 'villi'
adata_17025v.obs['smoking']= 'No'
adata_17025v.obs['cohort']= 'SP014'
adata_17025v.obs['library']= '10X v2'
adata_17025v.obs['procurement']= 'Graz'
adata_17025v.obs['placental_volume']= 'NA'
adata_17025v.obs['gestational_weeks']= 7 
adata_17025v.obs['gestational_days']= 54 
adata_17025v.obs['maternal_BMI']= 'NA' #info unavailable. 
adata_17025v.obs['maternal_age']= 'NA' #info unavailable. 


#Read the Cellbender corrected data for SP014 early: SEKR2-v. 
adata_sekr2v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_SEKR2-Villi_cellbender_filtered.h5')
adata_sekr2v.var_names_make_unique() 
adata_sekr2v.obs['donor_id']= 'Donor-SEKR2'
adata_sekr2v.obs['time']= 'early'
adata_sekr2v.obs['disease']= 'C'
adata_sekr2v.obs['tissue']= 'villi'
adata_sekr2v.obs['smoking']= 'No'
adata_sekr2v.obs['cohort']= 'SP014'
adata_sekr2v.obs['library']= '10X v2'
adata_sekr2v.obs['procurement']= 'Graz'
adata_sekr2v.obs['placental_volume']= 'NA'
adata_sekr2v.obs['gestational_weeks']= 7 
adata_sekr2v.obs['gestational_days']= 49
adata_sekr2v.obs['maternal_BMI']= 'NA' #info unavailable. 
adata_sekr2v.obs['maternal_age']= 'NA' #info unavailable. 

#Read the Cellbender corrected data for SP014 early: SOZE1-v. 
adata_soze1v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_SOZE1-Villi_cellbender_filtered.h5')
adata_soze1v.var_names_make_unique() 
adata_soze1v.obs['donor_id']= 'Donor-SOZE1'
adata_soze1v.obs['time']= 'early'
adata_soze1v.obs['disease']= 'C'
adata_soze1v.obs['tissue']= 'villi'
adata_soze1v.obs['smoking']= 'No'
adata_soze1v.obs['cohort']= 'SP014'
adata_soze1v.obs['library']= '10X v2'
adata_soze1v.obs['procurement']= 'Graz'
data_soze1v.obs['placental_volume']= 'NA'
adata_soze1v.obs['gestational_weeks']= 9 
adata_soze1v.obs['gestational_days']= 63 
adata_soze1v.obs['maternal_BMI']= 'NA' #info unavailable. 
adata_soze1v.obs['maternal_age']= 'NA' #info unavailable.  

#Rename the donor_id:
adata_soze1v.obs['donor_id']= 'Donor-SOZE1-villi'
adata_sekr2v.obs['donor_id']= 'Donor-SEKR2-villi'
adata_17025v.obs['donor_id']= 'Donor-17025-villi'

#Add few additional information to metadata:
#D11: 
adata_d11.obs['cohort']= 'SP082'
adata_d11.obs['library']= '10X v3'
adata_d11.obs['procurement']= 'Graz'
#D12: 
adata_d12.obs['cohort']= 'SP082'
adata_d12.obs['library']= '10X v3'
adata_d12.obs['procurement']= 'Graz'
#D34:
adata_d34.obs['cohort']= 'SP082'
adata_d34.obs['library']= '10X v3'
adata_d34.obs['procurement']= 'Graz'
#D39:
adata_d39.obs['cohort']= 'SP082'
adata_d39.obs['library']= '10X v3'
adata_d39.obs['procurement']= 'Graz'
#D41:
adata_d41.obs['cohort']= 'SP082'
adata_d41.obs['library']= '10X v3'
adata_d41.obs['procurement']= 'Graz'
#D37:
adata_d37.obs['cohort']= 'SP082'
adata_d37.obs['library']= '10X v3'
adata_d37.obs['procurement']= 'Graz'
#D38:
adata_d38.obs['cohort']= 'SP082'
adata_d38.obs['library']= '10X v3'
adata_d38.obs['procurement']= 'Graz'


#Merge all the early datasets into one object.
adata_early= adata_d11.concatenate(adata_d12, adata_d34, adata_d37, adata_d38, adata_d39, adata_d41, 
                                   adata_17025v, adata_sekr2v, adata_soze1v)


print(adata_early.obs['donor_id'].value_counts()) 



#Read SP014 late term control sample: Term-01. 
adata_term01 = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_Term-1_cellbender_filtered.h5')
adata_term01.var_names_make_unique() 
adata_term01.obs['donor_id']= 'Donor-Term01'
adata_term01.obs['time']= 'late_term'
adata_term01.obs['disease']= 'C'
adata_term01.obs['tissue']= 'villi'
adata_term01.obs['smoking']= 'No'
adata_term01.obs['cohort']= 'SP014'
adata_term01.obs['library']= '10X v2'
adata_term01.obs['procurement']= 'Graz'
adata_term01.obs['placental_volume']= 'NA'
adata_term01.obs['gestational_weeks']= 38 
adata_term01.obs['gestational_days']= 266 
adata_term01.obs['maternal_BMI']= 25.7 #info made available later (April'22)
adata_term01.obs['maternal_age']= 38  #info made available later (April'22)

#Read SP014 late term control sample: Term-02. 
adata_term02 = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_Term-2_cellbender_filtered.h5')
adata_term02.var_names_make_unique() 
adata_term02.obs['donor_id']= 'Donor-Term02'
adata_term02.obs['time']= 'late_term'
adata_term02.obs['disease']= 'C'
adata_term02.obs['tissue']= 'villi'
adata_term02.obs['smoking']= 'No'
adata_term02.obs['cohort']= 'SP014'
adata_term02.obs['library']= '10X v2'
adata_term02.obs['procurement']= 'Graz'
adata_term02.obs['placental_volume']= 'NA'
adata_term02.obs['gestational_weeks']= 38 
adata_term02.obs['gestational_days']= 267 
adata_term02.obs['maternal_BMI']= 28.9 #info made available later (April'22)
adata_term02.obs['maternal_age']= 32 #info made available later (April'22)

#Read SP014 late term control sample: Term-03. 
adata_term03 = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_Term-3_cellbender_filtered.h5')
adata_term03.var_names_make_unique() 
adata_term03.obs['donor_id']= 'Donor-Term03'
adata_term03.obs['time']= 'late_term'
adata_term03.obs['disease']= 'C'
adata_term03.obs['tissue']= 'villi'
adata_term03.obs['smoking']= 'No'
adata_term03.obs['cohort']= 'SP014'
adata_term03.obs['library']= '10X v2'
adata_term03.obs['procurement']= 'Graz'
adata_term03.obs['placental_volume']= 'NA'
adata_term03.obs['gestational_weeks']= 38 
adata_term03.obs['gestational_days']= 270 #266  
adata_term03.obs['maternal_BMI']= 20.2 #info made available later (April'22)
adata_term03.obs['maternal_age']= 33 #info made available later (April'22)



#Read SP014 late term control sample: 327-villi. 
adata_327v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_327-Placenta_cellbender_filtered.h5')
adata_327v.var_names_make_unique() 
adata_327v.obs['donor_id']= 'Donor-327-villi'
adata_327v.obs['time']= 'late_term'
adata_327v.obs['disease']= 'C'
adata_327v.obs['tissue']= 'villi'
adata_327v.obs['smoking']= 'No'
adata_327v.obs['cohort']= 'SP014'
adata_327v.obs['library']= '10X v2'
adata_327v.obs['procurement']= 'Oslo'
adata_327v.obs['placental_volume']= 'NA'
adata_327v.obs['gestational_weeks']= 39
adata_327v.obs['gestational_days']= 275
adata_327v.obs['maternal_BMI']= 22.68  
adata_327v.obs['maternal_age']= 35



#Read SP014 late term control sample: 328-villi. 
adata_328v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_328-Placenta_cellbender_filtered.h5')
adata_328v.var_names_make_unique() 
adata_328v.obs['donor_id']= 'Donor-328-villi'
adata_328v.obs['time']= 'late_term'
adata_328v.obs['disease']= 'C'
adata_328v.obs['tissue']= 'villi'
adata_328v.obs['smoking']= 'No'
adata_328v.obs['cohort']= 'SP014'
adata_328v.obs['library']= '10X v2'
adata_328v.obs['procurement']= 'Oslo'
adata_328v.obs['placental_volume']= 'NA'
adata_328v.obs['gestational_weeks']= 40 
adata_328v.obs['gestational_days']= 278
adata_328v.obs['maternal_BMI']= 25.78 
adata_328v.obs['maternal_age']= 38


#Read SP014 late term control sample: 372-villi. 
adata_372v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_372-Placenta_cellbender_filtered.h5')
adata_372v.var_names_make_unique() 
adata_372v.obs['donor_id']= 'Donor-372-villi'
adata_372v.obs['time']= 'late_term'
adata_372v.obs['disease']= 'C'
adata_372v.obs['tissue']= 'villi'
adata_372v.obs['smoking']= 'No'
adata_372v.obs['cohort']= 'SP014'
adata_372v.obs['library']= '10X v2'
adata_372v.obs['procurement']= 'Oslo'
adata_372v.obs['placental_volume']= 'NA'
adata_372v.obs['gestational_weeks']= 40 
adata_372v.obs['gestational_days']= 283
adata_372v.obs['maternal_BMI']= 23.05 
adata_372v.obs['maternal_age']= 37

#Concatenate all late term controls: 
adata_late= adata_term01.concatenate(adata_term02, adata_term03, adata_327v, adata_328v, adata_372v)
print(adata_late.obs['donor_id].value_counts()) 


#Merge the early & late into one object:
adata_control= adata_early.concatenate(adata_late)


#Read the SP014 PE data: 389 Villi 
adata_389v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_389-Placenta_cellbender_filtered.h5')
adata_389v.var_names_make_unique() 
adata_389v.obs['donor_id']= 'Donor-389-villi'
adata_389v.obs['time']= 'late_preterm'
adata_389v.obs['disease']= 'PE'
adata_389v.obs['tissue']= 'villi'
adata_389v.obs['smoking']= 'No'
adata_389v.obs['cohort']= 'SP014'
adata_389v.obs['library']= '10X v2'
adata_389v.obs['procurement']= 'Oslo'
adata_389v.obs['placental_volume']= 'NA' #placental wt= 220 
adata_389v.obs['gestational_weeks']= 27 
adata_389v.obs['gestational_days']= 192
adata_389v.obs['maternal_BMI']= 23.67 
adata_389v.obs['maternal_age']= 26

#Read the SP014 PE data: 419 villi 
adata_419v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_419-Placenta_cellbender_filtered.h5')
adata_419v.var_names_make_unique() 
adata_419v.obs['donor_id']= 'Donor-419-villi'
adata_419v.obs['time']= 'late_preterm'
adata_419v.obs['disease']= 'PE'
adata_419v.obs['tissue']= 'villi'
adata_419v.obs['smoking']= 'No'
adata_419v.obs['cohort']= 'SP014'
adata_419v.obs['library']= '10X v2'
adata_419v.obs['procurement']= 'Oslo'
adata_419v.obs['placental_volume']= 'NA' #placental wt= 280. 
adata_419v.obs['gestational_weeks']= 34 
adata_419v.obs['gestational_days']= 236
adata_419v.obs['maternal_BMI']= 21.88
adata_419v.obs['maternal_age']= 33


#Read the SP014 PE data: 557 villi (technical replicates: 557_1 & 557_2)
adata_557_1v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_577-Placenta_1_cellbender_filtered.h5')
adata_557_1v.var_names_make_unique() 
adata_557_1v.obs['donor_id']= 'Donor-557_1-villi'
adata_557_1v.obs['time']= 'late_preterm'
adata_557_1v.obs['disease']= 'PE'
adata_557_1v.obs['tissue']= 'villi'
adata_557_1v.obs['smoking']= 'No'
adata_557_1v.obs['cohort']= 'SP014'
adata_557_1v.obs['library']= '10X v2'
adata_557_1v.obs['procurement']= 'Oslo'
adata_557_1v.obs['placental_volume']= 'NA' 
adata_557_1v.obs['gestational_weeks']= 34 
adata_557_1v.obs['gestational_days']= 236 
adata_557_1v.obs['maternal_BMI']= 26.08
adata_557_1v.obs['maternal_age']= 35

#technical replicate: 557_2
adata_557_2v = sc.read_10x_h5('/data/analysis/preeclampsia_2019/analysis/images/SP014_cellbender_matrix/SP014_577-Placenta_2_cellbender_filtered.h5')
adata_557_2v.var_names_make_unique() 
adata_557_2v.obs['donor_id']= 'Donor-557_2-villi'
adata_557_2v.obs['time']= 'late_preterm'
adata_557_2v.obs['disease']= 'PE'
adata_557_2v.obs['tissue']= 'villi'
adata_557_2v.obs['smoking']= 'No'
adata_557_2v.obs['cohort']= 'SP014'
adata_557_2v.obs['library']= '10X v2'
adata_557_2v.obs['procurement']= 'Oslo'
adata_557_2v.obs['placental_volume']= 'NA' 
adata_557_2v.obs['gestational_weeks']= 34 
adata_557_2v.obs['gestational_days']= 236  
adata_557_2v.obs['maternal_BMI']= 26.08
adata_557_2v.obs['maternal_age']= 35

#Merge all early, late & PE to the control object:
adata= adata_control.concatenate(adata_389v, adata_419v, adata_557_1v, adata_557_2v)
print(adata.obs['donor_id'].value_counts()) 

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
plt.rcParams["figure.figsize"] = (15,15)

#Before filtering: 
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'donor_id', rotation= 90)

plt.rcParams["figure.figsize"] = (6,4)

#Before filtering: group by condition/time/library/procurement:
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'disease', rotation= 90)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'time', rotation= 90)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'library', rotation= 90)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0, groupby = 'procurement', rotation= 90)


#Filtering: A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. Here we will only consider cells with at least 200 detected genes and genes need to be expressed in at least 3 cells. Please note that those values are highly dependent on the library preparation method used.
adata_new= adata.copy()
#Perform all filtering on adata_new: 
sc.pp.filter_cells(adata_new, min_genes=200)
sc.pp.filter_genes(adata_new, min_cells=3)
print(adata_new.n_obs, adata_new.n_vars)


# Extremely high number of detected genes could indicate doublets. However, depending on the cell type composition in your sample, you may have cells with higher number of genes (and also higher counts) from one cell type. In this case, we will run doublet prediction further down, so we will skip this step now, but the code below is an example of how it can be run:
# Additionally, we can also see which genes contribute the most to such reads. We can for instance plot the percentage of counts per gene.
# P.S: Please don't do any technical doublet correction here. It might severely remove the SCT & EVT.
sc.pl.highest_expr_genes(adata, n_top=30) #MALAT1 is an ubiquitous gene. 

#filter for percent mito: keep 5% cut-off. (stick to this)
#The old data will be reanalyzed with 5% cut-off. 
adata_mito_filtered = adata_new[adata_new.obs['pct_counts_MT_genes'] < 5, :]

#After filtering: important QC plot. 
sc.pl.violin(adata_mito_filtered, ['pct_counts_MT_genes','pct_counts_Ribo_genes'], jitter=0, groupby = 'donor_id', rotation= 90)
sc.pl.violin(adata_mito_filtered, ['pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'disease', rotation= 90)
sc.pl.violin(adata_mito_filtered, ['pct_counts_MT_genes','pct_counts_Ribo_genes'],jitter=0, groupby = 'time', rotation= 90)

#After filtering: 
sc.pl.violin(adata_mito_filtered, ['n_genes_by_counts', 'total_counts', 'pct_counts_MT_genes','pct_counts_Ribo_genes'],
             jitter=0, groupby = 'donor_id', rotation= 90, save= '_placenta_all_QC_270122.pdf')
			 
sc.pl.violin(adata_mito_filtered, ['n_genes_by_counts', 'total_counts', 'pct_counts_MT_genes','pct_counts_Ribo_genes'], jitter=0, groupby = 'disease', rotation= 90, save= '_placenta_disease_QC_270122.pdf')


#XIST counts per cell: 
adata_mito_filtered.obs["XIST-counts"] = adata_mito_filtered.X[:,adata_mito_filtered.var_names.str.match('XIST')].toarray()
sc.pl.violin(adata_mito_filtered, ['XIST-counts'],jitter=0, groupby = 'donor_id', rotation= 90)


#Make scatter plots to see the relationship between total counts per cell vs #genes having one positive count per cell.
#Vienna recommendation: Don't blindly remove cells by "total counts" => might remove the tetraploid EVTs. 
sc.pl.scatter(adata_mito_filtered, x='total_counts', y='pct_counts_MT_genes')
sc.pl.scatter(adata_mito_filtered, x='total_counts', y='n_genes_by_counts')


#Save the filtered raw matrix
results= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_anndata/Placenta_cellbender_raw_matrix_270122.h5ad"
adata_mito_filtered.write(results)


#Cell-cycle scores: 
#Download from Github: https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
cell_cycle_genes = [x.strip() for x in open('/data/analysis/preeclampsia_2019/herse4olivia/regev_lab_cell_cycle_genes.txt')]
print(len(cell_cycle_genes))
#Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata_mito_filtered.var_names]
print(len(cell_cycle_genes))


#Before running cell cycle we have to normalize the data. In the scanpy object, the data slot will be overwritten with the normalized data. So first, save the raw data into the slot raw.
# Then run normalization, logarimize and scale the data.
adata_final= adata_mito_filtered.copy()
#Save the raw data safely in the 'raw' slot. 
adata_final.raw = adata_final

#Normalize to depth 10,000: 
sc.pp.normalize_per_cell(adata_final, counts_per_cell_after=1e4)
#take the logarithm:
sc.pp.log1p(adata_final)
#Scale data: scVI requires raw data, so scaling required. 

#Plot CC scores: 
sc.tl.score_genes_cell_cycle(adata_final, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pl.violin(adata_final, ['S_score', 'G2M_score'],jitter=0, groupby = 'donor_id', rotation=90)


#Make a count layer to securely store the raw data for subsequent scVI harmonization: 
adata_final.layers["counts"] = adata_final.raw.X.copy()

#Re-write the Cellbender processed, scanpy QC filtered & log-normalized object for SP014/SP082 object:  
latest= "/data/analysis/preeclampsia_2019/placenta_atlas_2022/villi_anndata/Placenta_cellbender_normalized_matrix_270122.h5ad"
adata_final.write(latest)




