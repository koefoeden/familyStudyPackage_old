library(ThodbergMisc)
library(Seurat)
library(Signac)
library(tidyverse)

library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

theme_set(theme_bw())

#### Files ####

# Aggregated sample
h5_fname <- "/Volumes/Groupdir/SUN-CBMR-SCOP_2021_0145/Output/220729_multiome_and_integration/aggregate/2022-08-23_aggregate/outs/filtered_feature_bc_matrix.h5"
frag_fname <- "/Volumes/Groupdir/SUN-CBMR-SCOP_2021_0145/Output/220729_multiome_and_integration/aggregate/2022-08-23_aggregate/outs/atac_fragments.tsv.gz"

#### Setup Signac/Seurat ####

# Gene counts
counts <- Read10X_h5(h5_fname)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
liver <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
)

# create ATAC assay and add it to the object
liver[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = frag_fname,
    annotation = annotation
)

# Create gene activtiy
GA <- GeneActivity(liver, assay = "ATAC")
liver[['PseudoRNA']] <- CreateAssayObject(counts = GA)

#### Calculate stats ####

DefaultAssay(liver) <- "ATAC"

# Nucleosome
liver <- NucleosomeSignal(liver)

# TSS enrichment
liver <- TSSEnrichment(liver)

# Counts in blacklist
liver$blacklist <- FractionCountsInRegion(liver, assay="ATAC", regions=blacklist_hg38_unified)

# Fragments in peaks
total_fragments <- CountFragments(frag_fname)
rownames(total_fragments) <- total_fragments$CB
liver$fragments <- total_fragments[colnames(liver), "frequency_count"]
liver <- FRiP(
    object = liver,
    assay = 'ATAC',
    total.fragments = 'fragments'
)

# Mitochrondia in genes
DefaultAssay(liver) <- "RNA"
liver[["percent.mt"]] <- PercentageFeatureSet(liver, pattern = "^MT-")

#### Add study ####

liver@meta.data$Biopsy <- colnames(liver) |>
    enframe(name=NULL) |>
    separate(value, into=c("Cell", "Study"), convert = TRUE) |>
    pull(Study)

#### Save ####

o_fname <- file.path(mount(), "emc/cbmr/users/nzl922/MultiomePrototype/liver_raw.qs")
qs::qsave(liver, file=o_fname, preset = "fast")

#### Tester @@@@

#VlnPlot(liver, features=colnames(liver@meta.data)[-1])

# subset(
#     x = liver,
#     subset = nFeature_RNA >= 1000 &
#         percent.mt < 1 &
#         nFeature_ATAC >= 1000 &
#         FRiP > 0.5 &
#         blacklist < 0.01 &
#         nucleosome_signal < 1.25 &
#         TSS.enrichment > 1
# )
