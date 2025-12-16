# Single-cell RNA-seq Analysis (PBMC 3k)
#### Parameter-robust scRNA-seq analysis using Seurat
(QC, PCA, clustering resolution stability validation)

## Project Overview
This project performs a single-cell RNA-seq analysis of PBMC 3k data with a focus on robust parameter selection.
Rather than relying on a single default setting, I systematically evaluated how QC thresholds, PCA dimensions, and clustering resolution affect downstream results, and selected biologically interpretable and stable conditions.

## Dataset
10x Genomics PBMC 3k
; Healthy donor peripheral blood mononuclear cells (scRNA-seq)
#### (pbmc3k/filtered_gene_bc_matrices/hg19)

## Tools & Environment
- R
- Seurat

## Analysis Workflow
Raw counts
 → Quality Control (nFeature, nCount, percent.mt)
 → Normalization
 → Highly Variable Gene selection
 → Scaling
 → PCA
 → Neighbor graph construction
 → Clustering
 → UMAP visualization
 → Marker identification & annotation
 #### (singlecell_best_set.R)

### A. QC validation (ver1_test_code.R)
#### Parameters test
- percent.mt: <8, <10, <15
- nFeature_RNA (upper bound): <2000, <2500, <3000
#### Conclusion
- mt < 10, nFeature < 2500 showed the most interpretable and stable UMAP structure
  (Overly strict or loose thresholds slightly degraded cluster clarity)

### B. PCA dimension selection (ver2_test_code.R)
#### Tested dimensions
PCs = 10 / 15 / 20
#### Conclusion
PCs = 10 (best balance between signal and interpretability)

### C. Clustering resolution Robustness (ver2_test_code.R)
#### Tested resolutions
res = 0.3 / 0.5 / 0.8
#### Conclusion
res = 0.5 provided optimal granularity without over-clustering

## Key Takeaways
- Single-cell analysis results can vary significantly with parameter choices
- Systematic validation is essential for robust and explainable results
- Marker consistency across conditions indicates true biological signal rather than parameter-driven artifacts

## Relevance
This project demonstrates:
- End-to-end scRNA-seq analysis using Seurat
- Practical QC and clustering validation strategies
- An analysis mindset suitable for NGS / single-cell data analysis services
