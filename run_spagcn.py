#!/bin/python

import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import time
import argparse
import anndata


def go_spagcn(adata):
    x_pixel = adata.obsm['spatial'][:, 0].tolist()
    y_pixel = adata.obsm['spatial'][:, 1].tolist()

    start_time = time.time()
    adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)

    p = 0.5
    l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    n_clusters = adata.obs['region'].unique().shape[0]
    #Set seed
    r_seed=t_seed=n_seed=100
    #Seaech for suitable resolution
    res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

    clf=spg.SpaGCN()
    clf.set_l(l)
    #Set seed
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    #Run
    clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
    y_pred, prob=clf.predict()
    adata.obs["pred"]= y_pred
    adata.obs["pred"]=adata.obs["pred"].astype('category')
    #Do cluster refinement(optional)
    #shape="hexagon" for Visium data, "square" for ST data.
    adj_2d=spg.calculate_adj_matrix(x=x_pixel,y=x_pixel, histology=False)
    refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="square")
    adata.obs["spaGCN"]=refined_pred
    adata.obs["spaGCN"]=adata.obs["spaGCN"].astype('category')

    end_time = time.time()
    processing_time = end_time - start_time

    adata.obs['spaGCN_time'] = processing_time
    return(adata)


def main():
    parser = argparse.ArgumentParser(description='Process some input parameters.')
    parser.add_argument('-i','--input_file', type=str, required=True, help='Path to the input file')
    parser.add_argument('-o','--output_file', type=str, required=True, help='Path to the output file')

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    
    print("loading data...")

    adata = anndata.read_h5ad(input_file)
    
    print("running spaGCN...")
    
    adata = go_spagcn(adata)
    
    print("outputing results...")
    
    adata.obs.to_csv(output_file)

if __name__ == '__main__':
    main()
