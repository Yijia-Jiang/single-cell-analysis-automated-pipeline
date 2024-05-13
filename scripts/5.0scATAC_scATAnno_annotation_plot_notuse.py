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
# parser.add_argument('--samplename', help="input name of the sample")
parser.add_argument('--input', help="input folder of count matrix from quickATAC")
# parser.add_argument('--ref', help="location of reference atlas")
# parser.add_argument('--atlas', help="indicate which reference atlas")
# parser.add_argument('--out', help="output directory")
parser.add_argument('--out_df', help="output directory")
#parser.add_argument('--uncertainty_threshold', help="uncertainty score threshold")
#parser.add_argument('--use_rep', help="use_rep")
#parser.add_argument('--distance_threshold', help="distance_threshold")
args = parser.parse_args()

sample_count_path=args.input
# reference_data_path=args.ref
# atlas = args.atlas
# out_h5ad=args.out
out_df=args.out_df

out_dir=os.path.dirname(out_h5ad)

output_name = args.samplename

print(sample_count_path) 
print(reference_data_path)
print(out_dir)

uncertainty_threshold = 0.5
distance_threshold = 95
reference_label_col = "celltypes"
use_rep = "X_spectral_harmony"

celltype_assignment_dir = os.path.join(out_dir, "celltype_assignment")
try:
    os.mkdir(celltype_assignment_dir)
except OSError as error:
    pass


query_annotated = sc.read_h5ad("sample_count_path")


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
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=KNN_pred_col, palette_dictionary=palette_dictionary, dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side" ), filename="1.query_celltype.png", title = "1. KNN-based celltypes")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=distance_pred_col,  dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="2.query_corrected_celltype.png", title = "2. Corrected celltypes")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=uncertainty_col, palette_dictionary=None, dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="2.query_corrected_uncertainty.png", title = "Uncertainty Score")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=cluster_col, palette_dictionary=None, dtype =None, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="3.query_leiden_number.png", title = "Leiden Clusters")
scATAnno_plotting.scATAnno_plotting_umap(query_annotated, hue=cluster_anno_col,  palette_dictionary=palette_dictionary, out_dir=os.path.join(celltype_assignment_dir,"query_assignment","Legend_side"), filename="3.query_cluster_annotation.png", title = "scATAnno Annotation")

    
