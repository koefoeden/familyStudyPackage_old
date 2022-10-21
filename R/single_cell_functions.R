#' @import Seurat
NULL

# Processing --------------------------------------------------------------

#' Create and save multimodal Seurat object
#'
#' @param filtered_h5_feature_bc_matrix_path
#' @param fragpath
#'
#' @return Object (and saves it as well)
#' @export
create_multimodal_object <- function(filtered_h5_feature_bc_matrix_path="data-raw/filtered_feature_bc_matrix.h5",
                                     frag_path = "data-raw/atac_fragments.tsv.gz",
                                     qs_save_path ="data/seur.qs") {

  counts <- Read10X_h5(filtered_h5_feature_bc_matrix_path)
  counts.rna <- counts %>% .$`Gene Expression`

  # get gene annotations for hg38
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"

  # create a Seurat object containing the RNA adata
  seur_obj <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
  )

  # create ATAC assay and add it to the object
  seur_obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = frag_path,
    annotation = annotation
  )

  qs::qsave(seur_obj, file = qs_save_path, preset = "fast")

  return(seur_obj)
}

#' Add Nucleosome signal, TSSEnrichment, Blacklist and peak counts, and mitochondiral genes
#'
#' @param seur_obj
#'
#' @return The seurat object with added stats
#' @export
add_stats_and_save <- function(seur_obj,
                               qs_save_path ="data/seur_stats.qs") {

  DefaultAssay(seur_obj) <- "ATAC"

  # Nucleosome
  seur_obj <- NucleosomeSignal(seur_obj)

  # TSS enrichment
  seur_obj <- TSSEnrichment(seur_obj)

  # Counts in blacklist
  seur_obj$blacklist <- FractionCountsInRegion(seur_obj, assay="ATAC", regions=blacklist_hg38_unified)

  # Fragments in peaks
  total_fragments <- CountFragments("data-raw/atac_fragments.tsv.gz")
  rownames(total_fragments) <- total_fragments$CB

  seur_obj$fragments <- total_fragments[colnames(seur_obj), "frequency_count"]

  seur_obj <- FRiP(
    object = seur_obj,
    assay = 'ATAC',
    total.fragments = 'fragments'
  )

  # Mitochrondia in genes
  DefaultAssay(seur_obj) <- "RNA"
  seur_obj[["percent.mt"]] <- PercentageFeatureSet(muscle, pattern = "^MT-")

  qs::qsave(seur_obj, file = qs_save_path, preset = "fast")
  return(seur_obj)
}

add_var_features_transform_dim_reducs_and_save <- function(seur_obj,
                                                           qs_save_path ="data/seur_stats_dimreducs.qs") {
  # Find Top Features RNA & ATAC
  seur_obj <- FindTopFeatures(seur_obj,
                              assay = "RNA", min.cutoff = 5)

  seur_obj <- FindTopFeatures(seur_obj,
                              assay="ATAC", min.cutoff = 5)

  # Transform RNA data and dimension reduce

  seur_obj <- SCTransform(seur_obj,
                          assay = "RNA",
                          vst.flavor="v2",
                          new.assay.name = "SCT")

  seur_obj <- RunPCA(seur_obj,
                     assay = "SCT",
                     reduction.name = "pca")

  # Dimension reduce ATAC data

  seur_obj <- RunTFIDF(seur_obj,
                       assay = "ATAC")

  seur_obj <- RunSVD(seur_obj,
                     assay = "ATAC", reduction.name = "lsi")

  qs::qsave(seur_obj, file = qs_save_path)

  return(seur_obj)
}

add_harmony_pca_and_lsi <- function(seur_obj) {
  DefaultAssay(seur_obj) <- 'ATAC'

  seur_obj <- RunHarmony(seur_obj, group.by.vars = "Biopsy",
                         reduction = "pca", dims.use = 1:30, reduction.save="harmony_pca", project.dim = FALSE,)

  seur_obj <- RunHarmony(seur_obj, group.by.vars = "Biopsy",
                         reduction = "lsi", dims.use = 1:30, reduction.save="harmony_lsi", project.dim = FALSE)

  return(seur_obj)
}

add_WNN_UMAP_clusters <- function(seur_obj,
                                  qs_save_path ="data/seur_stats_dimreducs_clusters.qs") {

  seur_obj <- FindMultiModalNeighbors(
    object = seur_obj,
    reduction.list = list("pca", "lsi"),
    dims.list = list(1:30, 2:30),
    verbose = TRUE)

  seur_obj <- RunUMAP(seur_obj, dims = 1:30, reduction = "pca", reduction.name = 'umap_rna')
  seur_obj <- RunUMAP(seur_obj, dims = 2:30, reduction = "lsi", reduction.name = 'umap_atac')
  seur_obj <- RunUMAP(seur_obj, nn.name = "weighted.nn", reduction.name = "umap_wnn")

  seur_obj <- FindClusters(seur_obj, graph.name = "wsnn", algorithm = 3, verbose = TRUE)

  qs::qsave(seur_obj, file = qs_save_path)

  return(seur_obj)
}
# QC ----------------------------------------------------------------------

plot_violin_plots_of_meta_data <- function(seur_object) {
  cell_data <- seur_object@meta.data %>%
    rownames_to_column(var="CB") %>%
    as_tibble() %>%
    select(-orig.ident)

  cell_data_long <- cell_data %>%
    pivot_longer(cols = 2:last_col())

  plot <- cell_data_long %>%
    ggplot(aes(x=factor(name), y=value)) +
    geom_violin() +
    facet_wrap(~name,
               scales="free",
               ncol=6) +
    geom_jitter(size=0.1, alpha=0.1) +
    ggtitle(str_glue("{ncol(seur_object)} cells"))

  return(plot)
}
