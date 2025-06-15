import pandas as pd
import anndata as ad
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap
import numpy as np
import scanpy as sc

# --- Configuration ---
# IMPORTANT: Replace with the actual path to your downloaded .h5ad file
H5AD_FILE_PATH = "/Users/venkatachalamsubramanianperiyasubbu/Downloads/8797c27c-6937-429e-818a-6f2bce18521a.h5ad" # <<<--- CHANGE THIS PATH!

# Name of the directory to store the output Parquet files
OUTPUT_DIR_NAME = "murine_glioblastoma_processed"

# PCA Configuration
N_PCA_COMPONENTS = 50 # Number of principal components to compute
APPLY_SCALING_BEFORE_PCA = True # Scale data before PCA

# UMAP Configuration
N_UMAP_COMPONENTS = 2 # Usually 2 or 3 for visualization
UMAP_N_NEIGHBORS = 15 # Adjust based on dataset size
UMAP_MIN_DIST = 0.1

# HVG Configuration
N_TOP_HVGS = 4000 # Number of highly variable genes to select. Adjust as needed.

# --- 1. Create Output Directory ---
os.makedirs(OUTPUT_DIR_NAME, exist_ok=True)
print(f"Created output directory: {OUTPUT_DIR_NAME}")

# --- 2. Load the H5AD file into an AnnData object ---
print("Loading AnnData object...")
try:
    adata = ad.read_h5ad(H5AD_FILE_PATH)
    print(f"Successfully loaded AnnData object from: {H5AD_FILE_PATH}")
    print(f"AnnData object shape: {adata.shape} (Cells x Genes)")
    print(f"Available layers: {list(adata.layers.keys())}")
    print(f"Available obsm keys: {list(adata.obsm.keys())}")
    print(f"Available varm keys: {list(adata.varm.keys())}")
    print(f"Available uns keys: {list(adata.uns.keys())}")

    # Ensure adata.X is dense for easier processing later if it's sparse
    if hasattr(adata.X, 'toarray'):
        adata.X = adata.X.toarray()
        print("Converted adata.X from sparse to dense array.")

except FileNotFoundError:
    print(f"Error: H5AD file not found at {H5AD_FILE_PATH}. Please check the path.")
    exit()
except Exception as e:
    print(f"An error occurred while loading the H5AD file: {e}")
    exit()

# --- 3. Basic Scanpy Preprocessing ---
print("\nStarting basic Scanpy preprocessing (QC, Normalization, Log-transformation)...")

# Basic filtering (adjust thresholds if needed after initial inspection)
print(f"  Initial cells: {adata.n_obs}, genes: {adata.n_vars}")
sc.pp.filter_cells(adata, min_genes=200) # Keep cells with at least 200 genes
sc.pp.filter_genes(adata, min_cells=3)   # Keep genes expressed in at least 3 cells
print(f"  After basic filtering: cells: {adata.n_obs}, genes: {adata.n_vars}")

# Normalize total counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)
# Log-transform the data
sc.pp.log1p(adata)
print("  Normalization and log-transformation complete.")


# --- 4. Save Core AnnData Components to Parquet ---
print("\nSaving core AnnData components to Parquet...")

# Save adata.X (Expression Matrix - now normalized and log-transformed)
expression_parquet_path = os.path.join(OUTPUT_DIR_NAME, "expression.parquet")
df_expression = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names) # adata.X should be dense now
df_expression.index.name = "cell_id"
df_expression.to_parquet(expression_parquet_path, index=True)
print(f"Saved expression data to: {expression_parquet_path}")

# Save adata.var (Gene Metadata)
gene_metadata_parquet_path = os.path.join(OUTPUT_DIR_NAME, "gene_metadata.parquet")
adata.var.index.name = "gene_id"
adata.var.to_parquet(gene_metadata_parquet_path, index=True)
print(f"Saved gene metadata to: {gene_metadata_parquet_path}")

# Save adata.obs (Cell Metadata)
cell_metadata_parquet_path = os.path.join(OUTPUT_DIR_NAME, "cell_metadata.parquet")
for col in adata.obs.select_dtypes(include=['category', 'object']).columns:
    adata.obs[col] = adata.obs[col].astype(str) # Convert categoricals/objects to string for parquet compatibility
adata.obs.index.name = "cell_id"
adata.obs.to_parquet(cell_metadata_parquet_path, index=True)
print(f"Saved cell metadata to: {cell_metadata_parquet_path}")


# --- 5. Perform PCA and save results ---
print(f"\nStarting PCA with {N_PCA_COMPONENTS} components...")

# Identify highly variable genes for dimensionality reduction
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=N_TOP_HVGS)
# Create a copy of adata subsetted to only highly variable genes for DR (reduces memory for PCA/UMAP)
adata_for_dr = adata[:, adata.var['highly_variable']].copy()
print(f"  Identified {adata_for_dr.shape[1]} highly variable genes for PCA.")

X_for_pca = adata_for_dr.X # Use the (now dense) subsetted HVG data

if APPLY_SCALING_BEFORE_PCA:
    print("  Scaling HVG data before PCA...")
    scaler = StandardScaler()
    X_for_pca = scaler.fit_transform(X_for_pca)
else:
    print("  Skipping scaling before PCA.")

pca = PCA(n_components=N_PCA_COMPONENTS, random_state=42)
pca_transformed_data = pca.fit_transform(X_for_pca)

pca_columns = [f"PC{i+1}" for i in range(N_PCA_COMPONENTS)]
df_pca = pd.DataFrame(pca_transformed_data, index=adata.obs_names, columns=pca_columns)
df_pca.index.name = "cell_id"

pca_parquet_path = os.path.join(OUTPUT_DIR_NAME, "pca_embeddings.parquet")
df_pca.to_parquet(pca_parquet_path, index=True)
print(f"Saved PCA embeddings to: {pca_parquet_path}")

df_explained_variance = pd.DataFrame({
    'PrincipalComponent': [f"PC{i+1}" for i in range(len(pca.explained_variance_ratio_))],
    'ExplainedVarianceRatio': pca.explained_variance_ratio_,
    'CumulativeExplainedVarianceRatio': np.cumsum(pca.explained_variance_ratio_)
})
explained_variance_parquet_path = os.path.join(OUTPUT_DIR_NAME, "pca_explained_variance.parquet")
df_explained_variance.to_parquet(explained_variance_parquet_path, index=False)
print(f"Saved PCA explained variance ratio to: {explained_variance_parquet_path}")


# --- 6. Perform UMAP and save results ---
print(f"\nStarting UMAP with {N_UMAP_COMPONENTS} components...")
X_for_umap = pca_transformed_data # UMAP typically runs on PCA results

reducer = umap.UMAP(n_components=N_UMAP_COMPONENTS,
                    n_neighbors=UMAP_N_NEIGHBORS,
                    min_dist=UMAP_MIN_DIST,
                    random_state=42,
                    transform_seed=42)

umap_embeddings = reducer.fit_transform(X_for_umap)

umap_columns = [f"UMAP{i+1}" for i in range(N_UMAP_COMPONENTS)]
df_umap = pd.DataFrame(umap_embeddings, index=adata.obs_names, columns=umap_columns)
df_umap.index.name = "cell_id"

umap_parquet_path = os.path.join(OUTPUT_DIR_NAME, "umap_embeddings.parquet")
df_umap.to_parquet(umap_parquet_path, index=True)
print(f"Saved UMAP embeddings to: {umap_parquet_path}")


# --- 7. Save Highly Variable Genes Metadata ---
if 'highly_variable' in adata.var.columns and adata.var['highly_variable'].any():
    df_hvg = adata.var[adata.var['highly_variable']].copy()
    hvg_metadata_parquet_path = os.path.join(OUTPUT_DIR_NAME, "highly_variable_gene_metadata.parquet")
    df_hvg.index.name = "gene_id"
    df_hvg.to_parquet(hvg_metadata_parquet_path, index=True)
    print(f"Saved highly variable gene metadata to: {hvg_metadata_parquet_path}")
else:
    print("No highly variable genes marked. Skipping saving HVG metadata.")


# --- 8. Save Basic Gene Statistics ---
print("\nCalculating basic gene statistics...")
df_gene_stats = pd.DataFrame(index=adata.var_names)
df_gene_stats.index.name = "gene_id"

# Use the normalized, log-transformed data for mean expression
df_gene_stats['mean_expression'] = np.asarray(adata.X).mean(axis=0)
# Count cells where expression > 0 (from current adata.X, which is log-transformed, so use >0)
df_gene_stats['n_cells_expressed'] = np.asarray(adata.X > 0).sum(axis=0)

gene_stats_parquet_path = os.path.join(OUTPUT_DIR_NAME, "gene_statistics.parquet")
df_gene_stats.to_parquet(gene_stats_parquet_path, index=True)
print(f"Saved gene statistics to: {gene_stats_parquet_path}")


print(f"\nAll essential Parquet files for Mouse Glioblastoma Atlas have been created in the '{OUTPUT_DIR_NAME}' directory.")
print("You can now use these files for your submission to the hackathon!")