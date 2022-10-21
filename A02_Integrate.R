library(ThodbergMisc)
library(Seurat)
library(Signac)
library(harmony)
library(tidyverse)

#### Load liver ####

# Read in data
#liver_fname <- file.path(ThodbergMisc::mount(), "emc/cbmr/users/nzl922/MultiomePrototype/aggregated_raw.rds")
liver_fname <- file.path(mount(), "emc/cbmr/users/nzl922/MultiomePrototype/liver_raw.qs")
liver0 <- qs::qread(liver_fname)

#### Filter ####

# Filter samples
liver <- subset(
    x = liver0,
    subset = nFeature_RNA >= 1000 &
        percent.mt < 1 &
        nFeature_ATAC >= 1000 &
        FRiP > 0.5 &
        blacklist < 0.01 &
        nucleosome_signal < 1.25 &
        TSS.enrichment > 1
)

# Filter genes
liver <- FindTopFeatures(liver, assay = "RNA", min.cutoff = 5)
liver <- FindTopFeatures(liver, assay="ATAC", min.cutoff = 5)

#### Reduce dimensions ####

liver <- SCTransform(liver, assay = "RNA", vst.flavor="v2", new.assay.name = "SCT")
liver <- RunPCA(liver, assay = "SCT", reduction.name = "pca")

liver <- RunTFIDF(liver, assay = "ATAC")
liver <- RunSVD(liver, assay = "ATAC", reduction.name = "lsi")

liver <- SCTransform(liver, assay = "PseudoRNA", vst.flavor="v2", new.assay.name = "PseudoSCT")
liver <- RunPCA(liver, assay = "PseudoSCT", reduction.name = 'pseudo_pca')

DefaultAssay(liver) <- 'ATAC'

#### Harmony ####

# Harmony
liver <- RunHarmony(liver, group.by.vars = "Biopsy", reduction = "pca", dims.use = 1:30, reduction.save="harmony_pca", project.dim = FALSE)
liver <- RunHarmony(liver, group.by.vars = "Biopsy", reduction = "lsi", dims.use = 1:30, reduction.save="harmony_lsi", project.dim = FALSE)
liver <- RunHarmony(liver, group.by.vars = "Biopsy", reduction = "pseudo_pca", dims.use = 1:30, reduction.save="harmony_pseudo", project.dim = FALSE)

#### WNN ####

# build a joint neighbor graph using both assays
# liver <- FindMultiModalNeighbors(
#     object = liver,
#     weighted.nn.name = "pre_wnn",
#     snn.graph.name = "pre_wsnn",
#     reduction.list = list("pca", "lsi"),
#     dims.list = list(1:30, 2:30),
#     verbose = TRUE)

liver <- FindMultiModalNeighbors(
    object = liver,
    #weighted.nn.name = "post_wnn",
    #snn.graph.name = "post_wsnn",
    reduction.list = list("harmony_pca", "harmony_lsi"),
    dims.list = list(1:30, 2:30),
    verbose = TRUE)

#### UMAP ####

# Pre
liver <- RunUMAP(liver, dims = 1:30, reduction = "pca", reduction.name = 'pre_rna')
liver <- RunUMAP(liver, dims = 2:30, reduction = "lsi", reduction.name = 'pre_atac')
liver <- RunUMAP(liver, dims = 1:30, reduction = "pseudo_pca", reduction.name = 'pre_pseudo')

# Post
liver <- RunUMAP(liver, dims = 1:30, reduction = "harmony_pca", reduction.name = 'umap_rna')
liver <- RunUMAP(liver, dims = 2:30, reduction = "harmony_lsi", reduction.name = 'umap_atac')
liver <- RunUMAP(liver, dims = 1:30, reduction = "harmony_pseudo", reduction.name = 'umap_pseudo')

# WNN
#liver <- RunUMAP(object = liver, nn.name = "pre_wnn", reduction.name = "pre_umap_wnn")
liver <- RunUMAP(object = liver,
                 nn.name = "weighted.nn",
                 reduction.name = "umap_wnn")

#### Cluster ####

liver <- FindClusters(liver, graph.name = "wsnn", algorithm = 3, verbose = TRUE)

#### Save ####

o_fname <- file.path(mount(), "emc/cbmr/users/nzl922/MultiomePrototype/liver_integrated.qs")
qs::qsave(liver, file=o_fname, preset = "fast")

#### Plot experiments ####

(DimPlot(liver, reduction = "pre_umap_rna", group.by = "Biopsy", shuffle = TRUE) + labs(title="RNA", x="UMAP1", y="UMAP2")) +
(DimPlot(liver, reduction = "pre_umap_wnn", group.by = "Biopsy", shuffle = TRUE) + labs(title="WNN", x="UMAP1", y="UMAP2")) +
(DimPlot(liver, reduction = "pre_umap_atac", group.by = "Biopsy", shuffle = TRUE) + labs(title="ATAC", x="UMAP1", y="UMAP2")) +
(DimPlot(liver, reduction = "post_umap_rna", group.by = "Biopsy", shuffle = TRUE) + labs(title="Harmony RNA", x="UMAP1", y="UMAP2")) +
(DimPlot(liver, reduction = "post_umap_wnn", group.by = "Biopsy", shuffle = TRUE) + labs(title="Harmony WNN", x="UMAP1", y="UMAP2")) +
(DimPlot(liver, reduction = "post_umap_atac", group.by = "Biopsy", shuffle = TRUE) + labs(title="Harmony ATAC", x="UMAP1", y="UMAP2")) +
    patchwork::plot_layout(guides = 'collect')

(DimPlot(liver, reduction = "pre_umap_rna", group.by = "seurat_clusters", label = TRUE) + labs(title="RNA", x="UMAP1", y="UMAP2")) +
    (DimPlot(liver, reduction = "pre_umap_wnn", group.by = "seurat_clusters", label = TRUE) + labs(title="WNN", x="UMAP1", y="UMAP2")) +
    (DimPlot(liver, reduction = "pre_umap_atac", group.by = "seurat_clusters", label = TRUE) + labs(title="ATAC", x="UMAP1", y="UMAP2")) +
    (DimPlot(liver, reduction = "post_umap_rna", group.by = "seurat_clusters", label = TRUE) + labs(title="Harmony RNA", x="UMAP1", y="UMAP2")) +
    (DimPlot(liver, reduction = "post_umap_wnn", group.by = "seurat_clusters", label = TRUE) + labs(title="Harmony WNN", x="UMAP1", y="UMAP2")) +
    (DimPlot(liver, reduction = "post_umap_atac", group.by = "seurat_clusters", label = TRUE) + labs(title="Harmony ATAC", x="UMAP1", y="UMAP2")) +
    patchwork::plot_layout(guides = 'collect')
#
# #### Label Transfer ####
#
# # Fetch atlases
# # macparland <- qs::qread("~/CBMR/emc/cbmr/users/nzl922/SymphonyReferences/Macparland/macparland_ref.qs")
# # macparland$save_uwot_path <- "~/CBMR/emc/cbmr/users/nzl922/SymphonyReferences/Macparland/macparland_uwot.uwot"
# #
# # # UMAP
# # query_macparland <- mapQuery(exp_query = liver$SCT@scale.data,
# #                   metadata_query = liver@meta.data,
# #                   ref_obj = macparland,
# #                   vars="Biopsy",
# #                   do_normalize = FALSE,
# #                   return_type="Seurat")
# #
# # # Cell type
# # query_macparland <- knnPredict.Seurat(query_macparland, macparland, 'Celltype')
# # query_macparland <- knnPredict.Seurat(query_macparland, macparland, 'Phase')
# #
# # liver$Celltype <- query_macparland$Celltype
# # liver$Phase <- query_macparland$Phase
#
# # Fetch atlases
# #atlas <- qs::qread("~/CBMR/emc/cbmr/users/nzl922/SymphonyReferences/LiverAtlas/liver_ref.qs")
# #atlas$save_uwot_path <- "~/CBMR/emc/cbmr/users/nzl922/SymphonyReferences/LiverAtlas/liver_uwot.uwot"
#
# atlas <- qs::qread("/Volumes/ELECOMSSD/Overflow/liver_ref.qs")
# atlas$save_uwot_path <- "/Volumes/ELECOMSSD/Overflow/liver_uwot.uwot"
#
#
# # UMAP
# query_atlas <- mapQuery(exp_query = liver$SCT@scale.data,
#                              metadata_query = liver@meta.data,
#                              ref_obj = atlas,
#                              vars="Biopsy",
#                              do_normalize = FALSE,
#                              return_type="Seurat")
#
# # Cell type
# query_atlas <- knnPredict.Seurat(query_atlas, atlas, 'cell_type')
# liver$Celltype <- droplevels(query_atlas$cell_type)
#
# #### Annotate and save ####
#
# # Annotations vs clusters
# cont_table <- table(liver$Celltype, liver$seurat_clusters)
#
# # Find dominant
# tmp <- apply(as.matrix(cont_table), 2, which.max)
# dominant_ct <- rownames(cont_table)[tmp]
#
# # Clean names
# dominant_ct <- factor(dominant_ct)
# dominant_ct <- dominant_ct |> fct_recode(LSECs="Endothelial Cells-LSECs",
#                           HSCs="Hepatic Stellate Cells and Fibroblasts",
#                           MAECs="Endothelial Cells-Macrovascular Arterial",
#                           MPs="Mononuclear Phagocytes",
#                           LymphoidCells="Lymphoid Cells") |>
#     as.character()
#
# new_labels  <- paste0(colnames(cont_table), "_", dominant_ct)
#
# liver$Cluster <- liver$seurat_clusters
# levels(liver$Cluster) <- new_labels
# Idents(liver) <- liver$Cluster
#
# DimPlot(liver, reduction = "post_umap_wnn", label = TRUE, repel=TRUE) + NoLegend()
#
# #### Save ####
#
# saveRDS(liver, file="/Volumes/ELECOMSSD/Overflow/aggregated_processed.rds")
#
#
# # Quick checks
# #DimPlot(liver, reduction = "post_umap_wnn", group.by = "Celltype", shuffle=TRUE, label = TRUE, repel=TRUE) + labs(title="Macparland", x="UMAP1", y="UMAP2")
# #DimPlot(liver, reduction = "post_umap_wnn", group.by = "Atlas", shuffle=TRUE, label = TRUE, repel=TRUE) + labs(title="Atlas", x="UMAP1", y="UMAP2")
# #DimPlot(liver, reduction = "post_umap_wnn", group.by = "Phase", shuffle=TRUE) + labs(title="Phase", x="UMAP1", y="UMAP2")
#
# # Sanity check
# # (DimPlot(liver, reduction = "post_umap_rna", group.by = "Celltype", label = TRUE) + labs(title="Harmony RNA", x="UMAP1", y="UMAP2")) +
# #     (DimPlot(liver, reduction = "post_umap_wnn", group.by = "Celltype", label = TRUE) + labs(title="Harmony WNN", x="UMAP1", y="UMAP2")) +
# #     (DimPlot(liver, reduction = "post_umap_atac", group.by = "Celltype", label = TRUE) + labs(title="Harmony ATAC", x="UMAP1", y="UMAP2")) +
# #     (DimPlot(liver, reduction = "post_umap_rna", group.by = "Atlas", shuffle=TRUE, label = TRUE, repel = TRUE) + labs(title="Harmony RNA", x="UMAP1", y="UMAP2")) +
# #     (DimPlot(liver, reduction = "post_umap_wnn", group.by = "Atlas", shuffle=TRUE, label = TRUE, repel = TRUE) + labs(title="Harmony WNN", x="UMAP1", y="UMAP2")) +
# #     (DimPlot(liver, reduction = "post_umap_atac", group.by = "Atlas", shuffle=TRUE, label = TRUE, repel = TRUE) + labs(title="Harmony ATAC", x="UMAP1", y="UMAP2")) +
# #     patchwork::plot_layout(guides = 'collect')
#
#
