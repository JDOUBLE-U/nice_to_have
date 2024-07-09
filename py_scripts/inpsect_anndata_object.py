#!/usr/bin/env python3

import anndata as ad

# Replace with the actual path to your .h5ad file
file_path = "/path/to/output_file.h5ad"

# Load the .h5ad file
adata = ad.read_h5ad(file_path)

# Print basic information about the loaded data
print(adata)

# Show the first few rows of the cell metadata (obs)
print(adata.obs.head())

# Show the first few rows of the gene metadata (var)
print(adata.var.head())

# Check the layers
print(adata.layers.keys())

# Example: Access spliced data
spliced_data = adata.layers['spliced']
print(spliced_data)