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

pbmc

# 4. 미토콘드리아 유전자 비율 계산
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data)

# 4-1. QC 시각화
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# 4-2. 필터 기준을 숫자로 설정
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10
)
pbmc

# 5. 정규화 (raw counts -> normalized data)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 6. 변동 유전자 선택 (고정보 유전자)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 6-1. 변동 유전자 확인 플롯
VariableFeaturePlot(pbmc)

# 7. 데이터 스케일링 (PCA 전에 보통 수행)
# percent.mt 같은 기술적 요인을 회귀로 제거할지 선택 가능
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# 8. PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc)

# 9. 이웃 그래프 생성 (PCA 공간에서 가까운 세포들 연결)
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# 10. 클러스터링 (해상도: 클러스터 개수 조절 노브)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 11. UMAP 임베딩 (시각화용 2D 좌표)
pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()

# 12. 클러스터별 마커 찾기 
markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

## 클러스터별 상위 마커 5개씩 확인
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

top_markers
print(top_markers, n = 50)


# 13. Labeling
new.cluster.ids <- c(
  "CD4 T (Naive)",
  "CD14+ Monocyte",
  "CD4 T (Memory)",
  "B cell",
  "CD8 T",
  "FCGR3A+ Monocyte",
  "NK cell",
  "Dendritic cell",
  "Platelet"
)

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

DimPlot(pbmc, reduction = "umap", label = TRUE) + NoLegend()



