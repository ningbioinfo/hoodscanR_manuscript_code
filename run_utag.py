#!/bin/bash

from utag import utag
import scanpy as sc
import anndata
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import argparse


def run_utag(adata, slide_key, resolution):
    start_time = time.time()
    utag_results = utag(
        adata,
        slide_key=slide_key,
        max_dist=15,
        normalization_mode='l1_norm',
        apply_clustering=True,
        clustering_method = 'leiden', 
        resolutions = [resolution]
    )
    end_time = time.time()
    processing_time = end_time - start_time
    utag_results.obs['utag_time'] = processing_time
    return(utag_results)

def main():
    parser = argparse.ArgumentParser(description='Process some input parameters.')
    parser.add_argument('-i','--input_file', type=str, required=True, help='Path to the input file')
    parser.add_argument('-o','--output_file', type=str, required=True, help='Path to the output file')
    parser.add_argument('-s','--slide_key', type=str, required=True, help='A string indicates the slide column name')
    parser.add_argument('-r','--resolution', type=float, required=True, help='A number, resolution for louvain')

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    slide_key = args.slide_key
    resolution = args.resolution
    
    print("loading data...")

    adata = anndata.read_h5ad(input_file)
    
    print("running Utag...")
    
    adata = run_utag(adata, slide_key, resolution)
    
    print("outputing results...")
    
    adata.obs.to_csv(output_file)

if __name__ == '__main__':
    main()






