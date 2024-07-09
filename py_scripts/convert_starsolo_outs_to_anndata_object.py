#!/usr/bin/env python3

import os
import pandas as pd
from scipy import io, sparse
import anndata as ad

def make_unique(var_names, join='-'):
    """
    Ensure unique variable names by appending a suffix.
    
    Parameters:
    var_names (pd.Index): The variable names to make unique.
    join (str): Delimiter to use for making variable names unique.
    
    Returns:
    pd.Index: Unique variable names.
    """
    seen = {}
    unique_var_names = []
    for name in var_names:
        if name in seen:
            seen[name] += 1
            unique_var_names.append(f"{name}{join}{seen[name]}")
        else:
            seen[name] = 0
            unique_var_names.append(name)
    return pd.Index(unique_var_names)

def starsolo_velocity_anndata(input_dir, output_file, delimiter='-'):
    """
    Convert STARsolo output into an AnnData object for RNA velocity analysis and save it to a file.

    Parameters:
    input_dir (str): Directory containing the STARsolo output files: barcodes.tsv, features.tsv, spliced.mtx, ambiguous.mtx, and unspliced.mtx
    output_file (str): Path to the output file where the AnnData object will be saved (.h5ad)
    delimiter (str): Delimiter to use for making variable names unique

    Returns:
    adata (AnnData): Annotated data matrix with spliced, ambiguous, and unspliced layers
    """
    # Load cell barcodes
    obs = pd.read_csv(os.path.join(input_dir, 'barcodes.tsv'), header=None, index_col=0)
    # Remove index column name to make it compliant with the AnnData format
    obs.index.name = None

    # Load gene features
    var = pd.read_csv(os.path.join(input_dir, "features.tsv"), sep='\t', header=None, names=['gene_ids', 'gene_names', 'feature_types'], index_col=1)
    var.index.name = None

    # Ensure unique variable names
    var.index = make_unique(var.index, join=delimiter)

    # Load matrix files
    spliced = sparse.csr_matrix(io.mmread(os.path.join(input_dir, "spliced.mtx")).T)
    ambiguous = sparse.csr_matrix(io.mmread(os.path.join(input_dir, "ambiguous.mtx")).T)
    unspliced = sparse.csr_matrix(io.mmread(os.path.join(input_dir, "unspliced.mtx")).T)

    # Create AnnData object
    adata = ad.AnnData(X=spliced, obs=obs, var=var, layers={'spliced': spliced, 'ambiguous': ambiguous, 'unspliced': unspliced})

    # Save the AnnData object to a file
    adata.write(output_file)

    return adata

# Usage :
adata = starsolo_velocity_anndata("/path/to/folder/raw", "/path/to/output.h5ad", delimiter='.')
