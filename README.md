
# snRNA-seq Data Preprocessing for Mouse Glioblastoma

This repository contains a Python script `processing.py` that processes a `.h5ad` single-nucleus RNA sequencing (snRNA-seq) dataset of mouse glioblastoma and converts it into optimized formats (`.parquet`) suitable for sharing and downstream analysis.

The **downloadable** dataset is available on Hugging Face:
[View on Hugging Face](https://huggingface.co/datasets/longevity-db/mouse-glioblastoma-snRNAseq)

## Dataset Summary

This dataset consists of single-nucleus RNA-sequencing (snRNA-seq) data of **mouse glioblastoma**. It is structured in the `.h5ad` format, commonly used in the AnnData ecosystem for storing large-scale omics data, particularly from [scanpy](https://scanpy.readthedocs.io/en/stable/).

## What the Script Does (`processing.py`)

The script performs a comprehensive preprocessing pipeline, converting the raw `.h5ad` file into more accessible and interoperable `.parquet` files for easier analysis and integration. It performs the following steps:

1.  **Setup Configuration**
    Defines the input path to the `.h5ad` file, sets the output directory name (`murine_glioblastoma_processed`), and configures parameters for PCA, UMAP, and Highly Variable Gene (HVG) selection.

2.  **Creates Output Directory**
    A directory named `murine_glioblastoma_processed` is created to hold all the output files, if it doesn't already exist.

3.  **Loads the `.h5ad` File**
    Reads the input `.h5ad` file into an AnnData object using the `anndata` library. It also converts the expression matrix (`adata.X`) to a dense array if it's sparse.

4.  **Basic Scanpy Preprocessing**
    Performs quality control filtering (e.g., minimum genes per cell, minimum cells per gene), normalizes total counts per cell to a target sum, and log-transforms the data.

5.  **Saves Core AnnData Components to Parquet**
    Converts and saves the following AnnData components as `.parquet` files:
    * `expression.parquet`: The normalized and log-transformed expression matrix (`adata.X`).
    * `gene_metadata.parquet`: Metadata about each gene (`adata.var`).
    * `cell_metadata.parquet`: Metadata about each cell (`adata.obs`), ensuring categorical/object columns are converted to strings for compatibility.

6.  **Perform PCA and save results**
    Identifies highly variable genes, scales the data (if configured), performs Principal Component Analysis (PCA) to reduce dimensionality, and saves the PCA embeddings (`pca_embeddings.parquet`) and explained variance ratios (`pca_explained_variance.parquet`).

7.  **Perform UMAP and save results**
    Applies Uniform Manifold Approximation and Projection (UMAP) on the PCA-transformed data to generate 2D embeddings for visualization, saving the results to `umap_embeddings.parquet`.

8.  **Save Highly Variable Genes Metadata**
    If highly variable genes are identified, their metadata is extracted from `adata.var` and saved to `highly_variable_gene_metadata.parquet`.

9.  **Save Basic Gene Statistics**
    Calculates the mean expression and the number of cells expressed for each gene (from the normalized, log-transformed data) and saves these statistics to `gene_statistics.parquet`.

## Output

The script saves the following `.parquet` files within the `murine_glioblastoma_processed` directory:

* `expression.parquet`: Gene expression matrix
* `gene_metadata.parquet`: Metadata about each gene
* `cell_metadata.parquet`: Metadata about each cell
* `pca_embeddings.parquet`: PCA embeddings (principal components for each cell)
* `pca_explained_variance.parquet`: Explained variance ratio for each principal component
* `umap_embeddings.parquet`: UMAP embeddings (2D or 3D coordinates for each cell)
* `highly_variable_gene_metadata.parquet`: Metadata for highly variable genes
* `gene_statistics.parquet`: Basic statistics for each gene (mean expression, number of cells expressed)

## Requirements

Install the required packages before running:

```bash
pip install pandas anndata scikit-learn umap-learn numpy scanpy pyarrow
```

## Usage

1.  **Download the `.h5ad` file:** Obtain the mouse glioblastoma `.h5ad` dataset (e.g., from the Hugging Face link provided above).
2.  **Update the script:** Replace the `H5AD_FILE_PATH` in the `processing.py` script with the actual path to your downloaded `.h5ad` file:

    ```python
    H5AD_FILE_PATH = "/path/to/your/mouse_glioblastoma_data.h5ad" # <<<--- CHANGE THIS PATH!
    ```
3.  **Run the script:** Execute the Python script from your terminal:

    ```bash
    python processing.py
    ```
4.  **Find the output:** The processed `.parquet` files will be saved in the `murine_glioblastoma_processed/` directory:

    ```
    murine_glioblastoma_processed/
    │
    ├── expression.parquet
    ├── gene_metadata.parquet
    ├── cell_metadata.parquet
    ├── pca_embeddings.parquet
    ├── pca_explained_variance.parquet
    ├── umap_embeddings.parquet
    ├── highly_variable_gene_metadata.parquet
    └── gene_statistics.parquet
    ```

## Notes

* Output files are in `.parquet` format, which is efficient for big data workflows and compatible with tools like Apache Spark, Pandas, and cloud platforms.

## References

- Source Dataset :
-     Original Data Source: CELLxGENE Discover Collection: "Glioblastoma from young and aged mice".
-     Direct .h5ad link: https://datasets.cellxgene.cziscience.com/8797c27c-6937-429e-818a-6f2bce18521a.h5ad
- Final Output: [https://huggingface.co/datasets/longevity-db/mouse-glioblastoma-snRNAseq](https://huggingface.co/datasets/longevity-db/mouse-glioblastoma-snRNAseq)
- File Format: [AnnData Documentation](https://anndata.readthedocs.io/en/latest/)
- Hugging Face Datasets: [https://huggingface.co/datasets](https://huggingface.co/datasets)

**Contributed by CellVPA Team Venkatachalam, Pooja, Albert**
