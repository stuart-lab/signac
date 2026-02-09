## code to prepare atac_small

library(Signac)
library(Seurat)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
setwd("vignette_data/pbmc_vignette/")
set.seed(1234)

ah <- AnnotationHub()
ensdb_v98 <- ah[["AH75011"]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

counts <- Read10X_h5(filename = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5")

gr_assay <- CreateGRangesAssay(
  counts = counts[1:100, 1:100],
  fragments = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz"
)

pbmc <- CreateSeuratObject(counts = gr_assay, assay = "peaks")

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

pbmc <- AddMotifs(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm[1:10])

Annotation(pbmc) <- annotations
ga <- GeneActivity(pbmc)

ga_subset <- ga[1:50, ]

pbmc[['RNA']] <- CreateAssayObject(counts = ga_subset)

annotation_subset <- annotations[seqnames(annotations) == 'chr1', ]
Annotation(pbmc) <- head(annotation_subset, 200)

pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc, features = rownames(pbmc))
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 1:10)
pbmc$cluster <- sample(c(1, 2), size = ncol(pbmc), replace = TRUE)

Fragments(pbmc) <- NULL
atac_small <- pbmc
usethis::use_data(atac_small, overwrite = TRUE)
