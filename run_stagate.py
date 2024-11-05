#!/bin/python

import STAGATE_pyG as STAGATE
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import anndata
import torch
import time
import argparse
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr

rpy2.robjects.numpy2ri.activate()
mclust = importr('mclust')

def go_stagate(adata, rad_cutoff, clust_num):
    start_time = time.time()
    STAGATE.Cal_Spatial_Net(adata, rad_cutoff=rad_cutoff)
    STAGATE.Stats_Spatial_Net(adata)
    adata = STAGATE.train_STAGATE(adata, device = 'cpu')
    sc.pp.neighbors(adata, use_rep='STAGATE')
    # clustering
    rmclust = robjects.r['Mclust']
    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(adata.obsm['STAGATE']), clust_num)
    mclust_res = np.array(res[-2])
    adata.obs['stagate'] = mclust_res.astype('int')
    end_time = time.time()
    processing_time = end_time - start_time
    adata.obs['stagate_time'] = processing_time
    return(adata)

def main():
    parser = argparse.ArgumentParser(description='Process some input parameters.')
    parser.add_argument('-i','--input_file', type=str, required=True, help='Path to the input file')
    parser.add_argument('-o','--output_file', type=str, required=True, help='Path to the output file')
    parser.add_argument('-d','--distance', type=int, required=True, help='A number describing the distance to be draw as neighborhood')
    #parser.add_argument('-r','--clust_num', type=float, required=True, help='A number, cluster number for mclust')

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    distance = args.distance
    
    print("loading data...")

    adata = anndata.read_h5ad(input_file)
    
    print("running stagate...")

    clust_num = len(adata.obs['region'].unique())
    
    adata = go_stagate(adata, distance, clust_num)
    
    print("outputing results...")
    
    adata.obs.to_csv(output_file)

if __name__ == '__main__':
    main()