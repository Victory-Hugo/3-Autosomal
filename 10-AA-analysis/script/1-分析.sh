#!/bin/bash
# 激活archetypal_analysis环境
# conde activate archetypal_analysis

# > archetypal-analysis -i INPUT_FILE -o OUTPUT_FILE -k K
# Arguments
# -i/--input_file: defines the input file / path. File must be in VCF, BED, PGEN or NPY format.
# -o/--output_file: defines the output file / path. File name does not need any extensions.
# -k/--n_archetypes: defines the number of archetypes.
# --tolerance: defines when to stop optimization.
# --max_iter: defines the maximum number of iterations.
# --random_state: defines the random seed number for initialization. No effect if "furthest_sum" is selected.
# -C/--constraint_coef: constraint coefficient to ensure that the summation of alfa's and beta's equals to 1. C is conisdered to be inverse of M^2 in the original paper.
# --initialize: defines the initialization method to guess initial archetypes.
# -dr/--dim_reduction: defines the dimensionality reduction technique to project the input data. Accepted=['PCA', 'MDS', 'UMAP', 'TSNE']. Default='PCA'.