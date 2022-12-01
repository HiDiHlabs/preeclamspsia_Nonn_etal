A repository for code used to analyze data and create figures presented in the following paper describing the single cell spatio-temporal landscape of 
placenta tissue in preeclampsia:
https://www.biorxiv.org/content/10.1101/2022.10.10.511539v1

Key contributors (in computational analysis):

1). CellRanger: Dr.Cornelius Fischer (Institute for Medical Systems Biology, Berlin)

2). CellBender: Dr. Cornelius Fischer & Olivia Debnath. 

3). snRNA seq analysis: Olivia Debnath (Berlin Institute of Health & Charité). The biological interpretation of cell type/state annotation was done 
   jointly with Daniela S. Valdes (MDC Berlin).

4). In-situ sequencing (ISS) analysis: Sebastian Tiesmeyer (Berlin Institute of Health & Charité). 

5). 10X visium spatial analysis: Kerim Ali Secener (MDC Berlin) with inputs from Olivia Debnath. 


Understanding the contents:
1). Placenta_celltyping: The scripts/notebooks dedicated towards data harmonization, batch-effect removal, cell type/state annotation, composition analysis, evaluation of integration & inferring trajectory of placenta villi.
2). Decidua_celltyping: The scripts/notebooks dedicated towards data harmonization, batch-effect removal, cell type/state annotation, composition analysis & evaluation of integration of maternal decidua.

3). RL_interaction_network: Scripts for inferring receptor-ligand interaction from 10X snRNA-seq data. The result files behind each figure generation were marked & explained.

4). 10X_visium_analysis: Scripts for analyzing 10X visium spatial data using placenta villi snRNA-seq as a reference for spot deconvolution.

5). Manuscript_analysis: Consists of notebooks/scripts to reproduce computational main & extended data figures of the manuscript.

Additionally, each subfolder has a "description" file to help the reader navigate through the contents and understand how each script/notebook function.
