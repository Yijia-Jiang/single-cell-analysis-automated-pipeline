from scipy.sparse import csr_matrix
import anndata as ad
import pandas as pd
import scipy.io
import os
import glob
import scanpy as sc
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--mtx', help="matrix.mtx")
parser.add_argument('--barcodes', help='barcodes.tsv')
parser.add_argument('--features', help='features.tsv')
parser.add_argument('--metadata', help='output')
parser.add_argument('--out', help='output')
# parser.add_argument('--features', help='features.tsv')
args = parser.parse_args()
mtx_file=args.mtx
cells_file=args.barcodes
features_file=args.features
meta_file=args.metadata
out_file=args.out
# # extend_range = "extend_2kb"
# promoter_peaks_file = "/mnt/cfce-rcsm/projects/development2/cfce_atac_annotation/BCC_GSE129785/software/QuickATAC/input/promoter_peaks/"+extend_range+ "/adult_peaks_promoters_" + extend_range + ".bed"
# output_name = "HealthyAult_simulated_9k_promoter_noharmony_" + extend_range
# out_dir = os.path.join("case_study", output_name)

# mtx_file = 'matrix.mtx'
# features_file = 'features.tsv'
# cells_file = 'barcodes.tsv'
# meta_file = 'metadata.csv'
# path='merged_sample_final'

def convert_mtx_to_h5ad(mtx_file,cells_file,features_file,meta_file,out_file):
  data = sc.read_mtx(os.path.join(mtx_file))
  data = data.T
  features = pd.read_csv(os.path.join(features_file), header=None)
  barcodes = pd.read_csv(os.path.join(cells_file),header=None)
  data.var_names = features[0]
  data.obs_names = barcodes[0]
  meta = pd.read_csv(os.path.join(meta_file), index_col=0)
  meta = meta.loc[data.obs_names,]
  data.obs = meta
  data.write(os.path.join(out_file))

convert_mtx_to_h5ad(mtx_file,cells_file,features_file,meta_file,out_file)
