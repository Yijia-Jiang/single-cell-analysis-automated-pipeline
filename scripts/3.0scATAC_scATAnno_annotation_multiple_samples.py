import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix
from scipy import sparse
import scipy.io
import os
import anndata as ad # Anndata version must > 0.8
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style='white')

import scATAnno
from scATAnno.SnapATAC2_spectral import *
from scATAnno.SnapATAC2_tools import *
from scATAnno.SnapATAC2_utils import *
from scATAnno import scATAnno_preprocess
from scATAnno import scATAnno_assignment
from scATAnno import scATAnno_integration
from scATAnno import scATAnno_plotting

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)

default_28 = scATAnno_plotting.get_palettes("default28")
default_102 = scATAnno_plotting.get_palettes("default102")

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--samplename', nargs="*", help="list of samples")
parser.add_argument('--input', nargs="*", help="list of input folders for count matrix from quickATAC")
parser.add_argument('--ref', help="data location of reference atlas")
parser.add_argument('--atlas', help="indicate which reference atlas")
parser.add_argument('--out', help="output directory")
parser.add_argument('--out_df', help="output annotation")
#parser.add_argument('--uncertainty_threshold', help="uncertainty score threshold")
#parser.add_argument('--use_rep', help="use_rep")
#parser.add_argument('--distance_threshold', help="distance_threshold")
args = parser.parse_args()
# sample_list =args.input
# print(sample_list)

# output_name_list = args.samplename
# print(output_name_list)


# sample_count_path=args.input
# input list of samples
output_name_list = args.samplename
sample_list =args.input
reference_data_path=args.ref
atlas = args.atlas
out_h5ad=args.out
out_df=args.out_df

out_dir=os.path.dirname(out_h5ad)

print(output_name_list )
print(sample_list)
# python xxxx.py --samplename 


# out_dir = os.path.join("case_study", output_name)


uncertainty_threshold = 0.5
distance_threshold = 95
reference_label_col = "celltypes"
use_rep = "X_spectral_harmony"
# atlas = "HealthyAdult"

# try:
#     os.mkdir("case_study")
# except OSError as error:
#     pass

# try:
#     os.mkdir(out_dir)
# except OSError as error:
#     pass

celltype_assignment_dir = os.path.join(out_dir, "celltype_assignment")
try:
    os.mkdir(celltype_assignment_dir)
except OSError as error:
    pass

import time 
import glob
start = time.time() 

# if atlas == "HealthyAdult":
# 	reference_data_path = "/mnt/cfce-rcsm/homes/yj976/projects/scATAnno/scATAnno-main/data/Reference/Healthy_Adult_reference_atlas.h5ad"
# elif atlas == "PBMC":
#     reference_data_path = "/mnt/cfce-rcsm/homes/yj976/projects/scATAnno/scATAnno-main/data/Reference/Healthy_Adult_reference_atlas.h5ad"
# elif atlas == "TIL":
#     reference_data_path = "/mnt/cfce-rcsm/homes/yj976/projects/scATAnno/scATAnno-main/data/Reference/Healthy_Adult_reference_atlas.h5ad"
# else:
#     print("Please specify reference atlas")

reference_data = scATAnno_preprocess.load_reference_data(reference_data_path)

data_list = []

for idx in range(len(sample_list)):
    query_data = scATAnno_preprocess.import_query_data(path = sample_list[idx],
                                    mtx_file = glob.glob(os.path.join(sample_list[idx], '*matrix.mtx.gz'))[0].split("/")[-1],
                                    cells_file = glob.glob(os.path.join(sample_list[idx],'*barcodes.tsv.gz'))[0].split("/")[-1],
                                    features_file = glob.glob(os.path.join(sample_list[idx],'*features.tsv.gz'))[0].split("/")[-1], 
                                    variable_prefix = output_name_list[idx], 
                                    celltype_col = "celltypes",
                                    add_metrics=False)
    data_list.append(query_data)

# for sample_count_path in sample_list:
#     query_data = scATAnno_preprocess.import_query_data(path = sample_count_path,
#                                     mtx_file = glob.glob(os.path.join(sample_count_path, '*matrix.mtx.gz'))[0].split("/")[-1],
#                                     cells_file = glob.glob(os.path.join(sample_count_path,'*barcodes.tsv.gz'))[0].split("/")[-1],
#                                     features_file = glob.glob(os.path.join(sample_count_path,'*features.tsv.gz'))[0].split("/")[-1], 
#                                     variable_prefix = output_name, 
#                                     celltype_col = "celltypes",
#                                     add_metrics=False)
#     data_list = query_data_list.append(data_list)

print("total number of query cells: {}".format(np.sum([i.obs.shape[0] for i in data_list])))

select_features(reference_data)
for anndata in data_list:
    select_features(anndata) 

datasets={}
datasets["Atlas"] = reference_data
for anndata in data_list:
    key = (anndata.obs["dataset"])[0]
    datasets[key] = anndata
datasets

end = time.time() 
print("Data Import time: " + str(end - start) +" seconds")

# Integration
# Integrate reference and query data
start = time.time() 
# Integrate reference and query data
integrated_adata = scATAnno_assignment.scATAnno_integrate_multiple(datasets, sample_size = 25000)

end = time.time() 
print("Dimension Reduction time: " + str(end - start) +" seconds")

integrated_adata.var_names_make_unique()
start = time.time() 
# Apply harmony to remove batch effects
integrated_adata_harmony = scATAnno_assignment.scATAnno_harmony(integrated_adata, batch_col = "dataset")
end = time.time() 
print("Harmony batch correction time: " + str(end - start) +" seconds")


start = time.time() 
# Plot UMAP using spectral embeddings
integrated_adata = scATAnno_assignment.scATAnno_umap(integrated_adata_harmony, out_dir, use_rep = use_rep, save = True)
end = time.time() 
print("UMAP time: " + str(end - start) +" seconds")



# integrated_adata = scATAnno_preprocess.load_integrated_adata(out_dir)
integrated_adata.var_names_make_unique()

start = time.time() 
reference = integrated_adata[integrated_adata.obs['dataset'] == "Atlas",:].copy()
reference.obs["celltypes"] =  scATAnno_assignment.curate_celltype_names(reference.obs["celltypes"], atlas = atlas)

# query Annotation
query = integrated_adata[integrated_adata.obs['dataset'] != "Atlas",:].copy()
# Perform KNN assignment
query_KNN = scATAnno_assignment.scATAnno_KNN_assign(reference, query, reference_label_col=reference_label_col, low_dim_col=use_rep)

# Perform weighted-distance based assignment
query_distance = scATAnno_assignment.scATAnno_distance_assign(reference, query_KNN, reference_label_col=reference_label_col, distance_threshold=distance_threshold, atlas=atlas, uncertainty_threshold=uncertainty_threshold, use_rep = use_rep)

# Perform cluster-level assignment
query_annotated = scATAnno_assignment.scATAnno_cluster_assign(query_distance, cluster_col = None, use_rep=use_rep)
end = time.time() 
print("Annotation time: " + str(end - start) +" seconds") # Annotation time: 61.26783514022827 seconds

# query_annotated.obs.to_csv(os.path.join(out_dir, "query_annotated.csv"))
query_annotated.obs.to_csv(os.path.join(out_df))
query_annotated.write(os.path.join(out_h5ad))
# query_annotated.write(os.path.join(out_dir, "2.query_annotated.h5ad"))

scATAnno_plotting.defaultPlotting_umap()
sc.pl.umap(integrated_adata, color="dataset",  show=True, title = "Integration")
ax = sc.pl.umap(integrated_adata, show=False)
# scATAnno_plotting.scATAnno_plotting_umap(integrated_adata[integrated_adata.obs.dataset!="Atlas"], hue="dataset",palette=['#228833'], title="Reference + Query", size=0.5, out_dir = os.path.join(celltype_assignment_dir, "projection"), filename="Integrated_Query.pdf", _ax=ax )
# scATAnno_plotting.scATAnno_plotting_umap(integrated_adata, hue="dataset",palette=['#228833','#cccccc'], size = 10, title="Reference + Query", out_dir = os.path.join(celltype_assignment_dir, "projection"), filename="Integrated_Query2.pdf")

scATAnno_plotting.scATAnno_plotting_umap(integrated_adata[integrated_adata.obs.dataset!="Atlas"], hue="dataset", title="Reference + Query", size=0.5, out_dir = os.path.join(celltype_assignment_dir, "projection"), filename="Integrated_Query.pdf", _ax=ax )
scATAnno_plotting.scATAnno_plotting_umap(integrated_adata, hue="dataset", size = 10, title="Reference + Query", out_dir = os.path.join(celltype_assignment_dir, "projection"), filename="Integrated_Query2.pdf")
scATAnno_plotting.scATAnno_plotting_umap(integrated_adata, hue="celltypes", size = 10, title="Reference + Query", out_dir = os.path.join(celltype_assignment_dir, "projection"), filename="Integrated_UMAP.pdf")

# query_annotated = sc.read_h5ad("merge_2.h5ad")
# out_dir = "./"

celltype_assignment_dir = os.path.join(out_dir, "celltype_assignment")


KNN_pred_col = '1.knn-based_celltype'
distance_pred_col = "2.corrected_celltype"
uncertainty_col='Uncertainty_Score'
cluster_col = 'leiden'
cluster_anno_col = 'cluster_annotation'
palette_dictionary = None

# Lengend on the side
try:
    os.mkdir(os.path.join(celltype_assignment_dir, "query_assignment","Legend_side"))
except OSError as error:
    pass
    
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=KNN_pred_col, palette_dictionary=palette_dictionary, dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side" ), filename="1.query_celltype.png", title = "1. KNN-based celltypes")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=distance_pred_col,  dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="2.query_corrected_celltype.png", title = "2. Corrected celltypes")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=uncertainty_col, palette_dictionary=None, dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="2.query_corrected_uncertainty.png", title = "Uncertainty Score")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=cluster_col, palette_dictionary=None, dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="3.query_leiden_number.png", title = "Leiden Clusters")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=cluster_anno_col,  palette_dictionary=palette_dictionary, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="3.query_cluster_annotation.png", title = "scATAnno Annotation")


scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=KNN_pred_col, palette_dictionary=palette_dictionary, dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side" ), filename="1.query_celltype.pdf", title = "1. KNN-based celltypes")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=distance_pred_col,  dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="2.query_corrected_celltype.pdf", title = "2. Corrected celltypes")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=uncertainty_col, palette_dictionary=None, dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="2.query_corrected_uncertainty.pdf", title = "Uncertainty Score")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=cluster_col, palette_dictionary=None, dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="3.query_leiden_number.pdf", title = "Leiden Clusters")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=cluster_anno_col,  palette_dictionary=palette_dictionary, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="3.query_cluster_annotation.pdf", title = "scATAnno Annotation")

# # Lengend in middle
# try:
#     os.mkdir(os.path.join(celltype_assignment_dir, "query_assignment","Legend_on"))
# except OSError as error:
#     pass

# sc.settings.autosave = False
# scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=KNN_pred_col, dtype ="categorical", palette_dictionary=palette_dictionary, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_on"), filename="1.celltype.pdf", title = "1. KNN-based celltypes")
# scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=distance_pred_col,  dtype ="categorical", palette_dictionary=palette_dictionary, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_on"), filename="2.corrected_celltype.pdf", title = "2. Corrected celltypes")
# scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=uncertainty_col, dtype ="numerical",  out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_on"), filename="2.corrected_uncertainty.pdf", title = "Uncertainty Score")
# scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=cluster_col, dtype ="categorical", out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_on"), filename="3.leiden_number.pdf", title = "Leiden Clusters")
# scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=cluster_anno_col, dtype ="categorical", palette_dictionary=palette_dictionary, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_on"), filename="3.cluster_annotation.pdf", title = "scATAnno Annotation")


# sc.settings.autosave = True
# sc.pl.umap(query_annotated, color = "cluster_annotation", title = "scATAnno Annotation", legend_loc="on data", save = "_on data_scATAnno_celltypes.pdf")
# sc.pl.umap(query_annotated, color ="Uncertainty_Score", save = "_scATAnno_uncertainty.pdf")

# sc.settings.autosave = False
# scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue="cluster_annotation", title="scATAnno Annotation",  out_dir = os.path.join(celltype_assignment_dir), filename="scATAnno_celltypes.pdf")
# scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue="Uncertainty_Score",  title = "scATAnno Uncertainty", out_dir = os.path.join(celltype_assignment_dir), filename="scATAnno_Uncertainty.pdf" )



