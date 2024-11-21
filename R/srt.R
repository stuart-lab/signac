library(patchwork)
library(tidyverse)
library(data.table)
library(Matrix.utils)
library(Seurat)
library(Signac)
library(tidyr)
library(dplyr)


pbmc.data <- Read10X(data.dir = " ")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


mito.data <- ReadMGATK(dir = " ")
# mito.data <- readRDS("~")
str(mito.data)
mito <- CreateAssayObject(counts = mito.data$counts)
mito <- subset(mito, cells = Cells(pbmc))
str(mito)
pbmc[["mito"]] <- mito
pbmc <- AddMetaData(pbmc, metadata = mito.data$depth[Cells(mito), ], col.name = "mtDNA_depth")
VlnPlot(pbmc, "mtDNA_depth", pt.size = 0.1) + scale_y_log10()
pbmc <- subset(pbmc, mtDNA_depth >= 0.1)
pbmc

DefaultAssay(pbmc) <- 'RNA'
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, label = TRUE) + NoLegend()

FeaturePlot(
  object = pbmc,
  features = c(" "),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 5
)
features<- c(" ")
dotplt = DotPlot(pbmc, features = features, cols = c('white','# '),dot.scale=4.5) + RotatedAxis()
dotplt = dotplt + theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 1))
print(dotplt)
pbmc <- RenameIdents(
  object = pbmc,
  '0' = ' ',
  '1' = ' '
)
FeatureScatter(pbmc, "mtDNA_depth", "nucleosome_signal") + ggtitle("") + scale_x_log10()

variable.sites <- IdentifyVariants(pbmc, assay = "mito", refallele = mito.data$refallele)
# variable.sites <- fread(cmd = paste("zcat", "~"))
VariantPlot(variants = variable.sites)
# proper method to draw plot !!!!!!
high.conf <- variable.sites[variable.sites$n_cells_conf_detected >= 3, ]
high.conf$pos <- high.conf$vmr > 0.01 & high.conf$strand_correlation > 0.65
ggplot(data = high.conf, mapping = aes_string(x = "strand_correlation", y = "vmr", color = "pos")) + 
    geom_point() + labs(x = "Strand concordance", y = "Variance-mean ratio") + 
    geom_vline(xintercept = 0.65, color = "black", linetype = 2) + 
    geom_hline(yintercept = 0.01, color = "black", linetype = 2) + 
    scale_color_manual(values = c("black", "firebrick")) + 
    scale_y_continuous(labels = comma) + theme_classic() + theme(legend.position = "none")

high.conf <- subset(
  variable.sites, subset = n_cells_conf_detected >= 3 &  
    strand_correlation >= 0.65 &
    vmr > 0.01
)
high.conf[,c(1,2,5)]
pbmc <- AlleleFreq(
  object = pbmc,
  variants = high.conf$variant,
  assay = "mito"
)
pbmc[["alleles"]]
DefaultAssay(pbmc) <- "alleles"
alleles.view <- c("14384G>A", "215A>G", "385A>G")
FeaturePlot(
  object = pbmc,
  features = alleles.view,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 3
) & NoLegend()
DoHeatmap(pbmc, features = rownames(pbmc), slot = "data", disp.max = 1) +
  scale_fill_viridis_c()
DimPlot(pbmc, label = TRUE) + NoLegend()
