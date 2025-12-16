#############################
#                           #
#   resolution, dims test   #
#                           #
#############################

setwd("D:/project_04_singlecell")

library(Seurat)
library(dplyr)

# 1. 데이터 다운로드 및 압축해제
dir.create("pbmc3k")
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz", exdir = "pbmc3k")

# 2. 데이터 로딩
pbmc.data <- Read10X(data.dir = "pbmc3k/filtered_gene_bc_matrices/hg19")

# 3. Seurat object 생성
pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "PBMC3K",
  min.cells = 3,
  min.features = 200
)

# 4. 기본 확인
pbmc

################################################################################

# 1. 미토콘드리아 유전자 비율 계산
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 1-1. 계산 확인
head(pbmc@meta.data)

# 1-2. QC 시각화
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# 1-3. 필터 기준을 숫자로 설정
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10
)
pbmc

################################################################################

# 1) 정규화 (raw counts -> normalized data)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 2) 변동 유전자 선택 (고정보 유전자)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

## 변동 유전자 확인 플롯(HVG)
VariableFeaturePlot(pbmc)

# 3) 데이터 스케일링 (PCA 전에 보통 수행)
# percent.mt 같은 기술적 요인을 회귀로 제거할지 선택 가능
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# 4) PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# PCA 결과 확인
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc)

pbmc_base <- pbmc

################################################################################

# 1. 데이터 다운로드 및 압축해제
dir.create("pbmc3k")
untar("pbmc3k_filtered_gene_bc_matrices.tar.gz", exdir = "pbmc3k")

# 2. 데이터 로딩
pbmc.data <- Read10X(data.dir = "pbmc3k/filtered_gene_bc_matrices/hg19")

# 3. Seurat object 생성
pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "PBMC3K",
  min.cells = 3,
  min.features = 200
)

# 4. 기본 확인
pbmc

################################################################################

# 1. 미토콘드리아 유전자 비율 계산
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 1-1. 계산 확인
head(pbmc@meta.data)

# 1-2. QC 시각화
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# 1-3. 필터 기준을 숫자로 설정
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10
)
pbmc

################################################################################

# 1) 정규화 (raw counts -> normalized data)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 2) 변동 유전자 선택 (고정보 유전자)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

## 변동 유전자 확인 플롯(HVG)
VariableFeaturePlot(pbmc)

# 3) 데이터 스케일링 (PCA 전에 보통 수행)
# percent.mt 같은 기술적 요인을 회귀로 제거할지 선택 가능
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# 4) PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# PCA 결과 확인
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc)

################################################################################
pbmc_base <- pbmc

library(ggplot2)

run_cluster_umap <- function(obj, dims_use, res = 0.5){
  obj <- FindNeighbors(obj, dims = 1:dims_use)
  obj <- FindClusters(obj, resolution = res)
  obj <- RunUMAP(obj, dims = 1:dims_use)
  return(obj)
}

# dims 비교 (res = 0.5 고정)
pbmc_d10 <- run_cluster_umap(pbmc_base, dims_use = 10, res = 0.5)
pbmc_d15 <- run_cluster_umap(pbmc_base, dims_use = 15, res = 0.5)
pbmc_d20 <- run_cluster_umap(pbmc_base, dims_use = 20, res = 0.5)

cat("dims10 clusters:", length(levels(Idents(pbmc_d10))), "\n")
cat("dims15 clusters:", length(levels(Idents(pbmc_d15))), "\n")
cat("dims20 clusters:", length(levels(Idents(pbmc_d20))), "\n")

p10 <- DimPlot(pbmc_d10, reduction = "umap", label = TRUE) + ggtitle("dims=10") + NoLegend()
p15 <- DimPlot(pbmc_d15, reduction = "umap", label = TRUE) + ggtitle("dims=15") + NoLegend()
p20 <- DimPlot(pbmc_d20, reduction = "umap", label = TRUE) + ggtitle("dims=20") + NoLegend()

library("patchwork")
(p10 | p15 | p20)

################################################################################

# resolution 비교 (dims = 10 고정)
pbmc_r03 <- run_cluster_umap(pbmc_base, dims_use = 10, res = 0.3)
pbmc_r05 <- run_cluster_umap(pbmc_base, dims_use = 10, res = 0.5)
pbmc_r08 <- run_cluster_umap(pbmc_base, dims_use = 10, res = 0.8)

cat("res0.3 clusters:", length(levels(Idents(pbmc_r03))), "\n")
cat("res0.5 clusters:", length(levels(Idents(pbmc_r05))), "\n")
cat("res0.8 clusters:", length(levels(Idents(pbmc_r08))), "\n")

p0.3 <- DimPlot(pbmc_r03, reduction="umap", label=TRUE) + ggtitle("dims=10, res=0.3") + NoLegend()
p0.5 <- DimPlot(pbmc_r05, reduction="umap", label=TRUE) + ggtitle("dims=10, res=0.5") + NoLegend()
p0.8 <- DimPlot(pbmc_r08, reduction="umap", label=TRUE) + ggtitle("dims=10, res=0.8") + NoLegend()

(p0.3 | p0.5 | p0.8)

get_top3 <- function(obj){
  m <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
  m %>% group_by(cluster) %>% slice_max(avg_log2FC, n=3)
}

top3_r03 <- get_top3(pbmc_r03)
top3_r05 <- get_top3(pbmc_r05)
top3_r08 <- get_top3(pbmc_r08)

################################################################################



