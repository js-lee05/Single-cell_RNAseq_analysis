#############################
#                           #
#        QC 기준 test       #
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

################################################################################

# A. percent.mt 임계값 비교
pbmc_raw <- pbmc

qc_grid <- list(
  mt8  = list(min_feat=200, max_feat=2500, mt=8),
  mt10 = list(min_feat=200, max_feat=2500, mt=10),
  mt15 = list(min_feat=200, max_feat=2500, mt=15)
)

run_qc <- function(obj, min_feat, max_feat, mt){
  subset(obj, subset = nFeature_RNA > min_feat & nFeature_RNA < max_feat & percent.mt < mt)
}

for (nm in names(qc_grid)) {
  p <- qc_grid[[nm]]
  tmp <- run_qc(pbmc_raw, p$min_feat, p$max_feat, p$mt)
  cat(nm, ": cells=", ncol(tmp), "features=", nrow(tmp), "\n")
}

################################################################################

# B. doublet 의심 상한선(max_feat) 비교
qc_grid2 <- list(
  max2000 = list(min_feat=200, max_feat=2000, mt=10),
  max2500 = list(min_feat=200, max_feat=2500, mt=10),
  max3000 = list(min_feat=200, max_feat=3000, mt=10)
)

for (nm in names(qc_grid2)) {
  p <- qc_grid2[[nm]]
  tmp <- run_qc(pbmc_raw, p$min_feat, p$max_feat, p$mt)
  cat(nm, ": cells=", ncol(tmp), "\n")
}

################################################################################

# 각 QC 조건에서 UMAP 구조 비교
pbmc_qc_base <- pbmc_raw 

run_full <- function(obj, min_feat=200, max_feat=2500, mt=10, dims_use=10, res=0.5){
  obj <- subset(obj, subset = nFeature_RNA > min_feat & nFeature_RNA < max_feat & percent.mt < mt)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method="vst", nfeatures=2000)
  obj <- ScaleData(obj, features=rownames(obj))
  obj <- RunPCA(obj, features=VariableFeatures(obj))
  obj <- FindNeighbors(obj, dims=1:dims_use)
  obj <- FindClusters(obj, resolution=res)
  obj <- RunUMAP(obj, dims=1:dims_use)
  return(obj)
}

# A) mt cutoff 비교
pbmc_mt8  <- run_full(pbmc_qc_base, mt=8,  max_feat=2500)
pbmc_mt10 <- run_full(pbmc_qc_base, mt=10, max_feat=2500)
pbmc_mt15 <- run_full(pbmc_qc_base, mt=15, max_feat=2500)

# B) max nFeature 비교
pbmc_max2000 <- run_full(pbmc_qc_base, mt=10, max_feat=2000)
pbmc_max2500 <- run_full(pbmc_qc_base, mt=10, max_feat=2500)
pbmc_max3000 <- run_full(pbmc_qc_base, mt=10, max_feat=3000)


library(patchwork)

(DimPlot(pbmc_mt8,  label=TRUE) + ggtitle("mt<8")  + NoLegend()) |
  (DimPlot(pbmc_mt10, label=TRUE) + ggtitle("mt<10") + NoLegend()) |
  (DimPlot(pbmc_mt15, label=TRUE) + ggtitle("mt<15") + NoLegend())

(DimPlot(pbmc_max2000, label=TRUE) + ggtitle("nFeature<2000") + NoLegend()) |
  (DimPlot(pbmc_max2500, label=TRUE) + ggtitle("nFeature<2500") + NoLegend()) |
  (DimPlot(pbmc_max3000, label=TRUE) + ggtitle("nFeature<3000") + NoLegend())


# 마커 안정성 비교
get_top3 <- function(obj){
  m <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
  m %>% group_by(cluster) %>% slice_max(avg_log2FC, n=3)
}

top3_mt8  <- get_top3(pbmc_mt8)
top3_mt10 <- get_top3(pbmc_mt10)
top3_mt15 <- get_top3(pbmc_mt15)

top3_max2000 <- get_top3(pbmc_max2000)
top3_max2500 <- get_top3(pbmc_max2500)
top3_max3000 <- get_top3(pbmc_max3000)
