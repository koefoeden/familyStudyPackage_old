
# *** CORRECTED get_MD_data_for sample instead of get_MA_data bug ***
#' Get MD data for all samples in long tibble
#'
#' @param object The DGElist object with transcript/gene counts and sample information
#' encoded in the y$sample object
#' @return A long tibble with three columns: sample label, mean and diff values
#' @export
#' @examples
#' get_MD_data_for_all_samples()

get_MD_data_for_all_samples_updated <- function(object) {
  sample_indices_w_names <- colnames(object) %>%
    purrr::set_names()

  MA_plot_data_long <- purrr::map(sample_indices_w_names,
                                  ~cbmr::get_MD_data_for_sample(object = object,
                                                                column =.)) %>%
    dplyr::bind_rows(.id="Sample") %>%
    dplyr::mutate(Sample=factor(Sample,
                                levels=gtools::mixedsort(colnames(object))))
}

# *** ADDED: Now with added cols option ****
#' Plot MD-figures for all samples in DGElist.
#'
#' @param object The DGElist object with transcript/gene counts and sample information
#' encoded in the y$sample object
#' @param samples The samples to plot as a character vector
#' @param ncol Number of columns in the combined plot as integer vector
#' @return A facetted ggplot
#' @export
#' @examples
#' ggplot_MD()
ggplot_MD_updated <-function(object, samples="all", ncol=3) {
  MD_data <- get_MD_data_for_all_samples_updated(object)

  if (!identical(samples, "all")) {
    MD_data <- MD_data %>%
      dplyr::filter(Sample %in% as.character(samples))
  }

  MD_data %>%
    ggplot(aes(x=Mean, y=Diff, col=Diff>0)) +
    geom_point(size=1, alpha=0.5) +
    facet_wrap(~Sample,ncol = ncol) +
    xlab("Average log CPM (this sample and others)") +
    ylab("log-ratio (this sample vs others)") +
    theme(legend.position = "none")
}

n_sig_genes_pr_contrast <- function(results_DF_list) {
  # *** NEW FUNCTION****
  #' Get number of significant genes/GO-terms per contrast
  #' @param results_DF_list An list of either DE genes or GO-terms results,
  #' with an element for each contrast.
  #' @return A facetted ggplot
  #' @export
  #' @examples
  #' ggplot_MD()
  purrr::map(results_DF_list,
             ~.x %>%
               filter(if_any(any_of(c("FDR", "adj.P.Val")),
                                   ~.x < 0.05)) %>%
               nrow()) %>%
    unlist()
}

# **** CHANGED: Added documentation & export, and IS NOW INTERACTIVE  ****
#' Make a volcano plot using ggplot
#'
#' @param df gene results from EdgeR or voom
#' @param genes character vector, "all" for plotting all genes,
#' otherwise only plot vector of genes
#' @return The ggplot object
#' @export
ggplot_volcano_updated <- function(df, genes="all") {
  if(!identical(genes,"all")) {
    df <- df %>% filter(Name %in% genes)
  }

  allLogFc <- df %>% pull("logFC")
  minLogFc <- allLogFc %>% min()
  maxLogFc <- allLogFc %>% max()

  minPval <- df %>%
    select(any_of(c("PValue", "P.Value"))) %>%
    pull() %>%
    log10() %>%
    `*`(-1) %>%
    max()

  maxPval <- 0

  df_formatted <- df %>%
    mutate(across(any_of(c("PValue", "P.Value")),
                  .fns = ~-log10(.x),
                  .names="log10Pval")) %>%
    mutate(across(any_of(c("FDR","adj.P.Val")),
                  .fns = ~.x < 0.05,
                  .names="Significant_after_adjustment")) %>%
    mutate(gene_is_interesting = Name %in% get_interesting_genes())

  plot <- ggplot(df_formatted, aes_string(x = "logFC",
                                          y = "log10Pval",
                                          colour = "Significant_after_adjustment",
                                          text="Name",
                                          alpha="gene_is_interesting")) +
    geom_point(size = 0.5) +
    scale_color_manual(values = c(`TRUE` = "red",
                                  `FALSE` = "black"),
                       name = "Significance",
                       labels = c(`TRUE` = "adj.P.Val < 0.05",
                                  `FALSE` = "adj.P.Val ≥ 0.05"),
                       breaks = c("TRUE", "FALSE")) +
    scale_alpha_manual(values = c(`TRUE` = 1,
                                  `FALSE` = 0.5)) +
    ylab("-log10 P-value") +
    xlab("Log2 Fold Change") +
    scale_y_continuous(expand = expansion(c(0,0.5))) +
    geom_hline(yintercept = -log10(0.05), lty="dashed", col="grey") +
    geom_vline(xintercept = 0, lty="dashed", col="grey")

  interactive_plot <- plot %>% ggplotly()
  return(interactive_plot)
}

# *** CHANGED: NOW EXPECTS BOTH GENE NAME AND ENSEMBL_ID***
#' Return a nicely formatted DT::datatable object of the DE results
#'
#' @param df The dataframe with the DE results
#' @param type Whether the results is edgeR or voom
#'
#' @return A DT::datatable
#' @export
#'
#' @examples
#' get_DE_datatable()
get_DE_datatable_updated <- function(df) {
  df %>%
    select(any_of(c("Name", "ENSEMBL_ID", "description", "logFC", "PValue","FDR")))  %>%
    mutate(across(any_of(c("logFC", "PValue","FDR")),
                  ~signif(.x,2))) %>%
    DT::datatable(extensions = 'Buttons',
                  filter="top",
                  rownames = FALSE,
                  options = list(dom = 'Bfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'colvis'),
                                 pageLength = 10),
                  height = "100%", width = "100%")
}

# *** CORRECTED: Can now take multiple color cols...***
#' Plot samples in two-dimensional space using MDS
#'
#' @param y The DGE list object
#' @param dims The dimensions to plot as a vector
#' @param color_by The vector to colour the plots by
#' @return Returns a list of ggplots coloured according to color_by
#' @export
ggplot_mds_repel_updated <- function (y, dims, color_by) {

  plotMDS_obj <- edgeR::plotMDS.DGEList(y, dim.plot = dims,
                                        plot = FALSE)
  x_y_data <- plotMDS_obj$eigen.vectors[, dims] %>% as.data.frame()
  x_y_data <- cbind(x_y_data, y$samples)

  var_explained_per_dim <- plotMDS_obj$var.explained[dims] %>%
    signif(2) %>% `*`(100)

  axis_labels <- plotMDS_obj$axislabel


  plots <- map(color_by,
               ~x_y_data %>%
                 ggplot(aes_string(x = colnames(x_y_data)[1],
                                   y = colnames(x_y_data)[2],
                                   text="Salmon.ID.edited",
                                   color = .x)) +
                 ggplot2::geom_point() +
                 xlab(str_glue("{axis_labels}. {dims[1]} ({var_explained_per_dim[1]} % var. explained)")) +
                 ylab(str_glue("{axis_labels}. {dims[2]} ({var_explained_per_dim[2]} % var. explained)")) +
                 ggtitle(str_glue("MDS-plot colored by {.x}. Dimensions: {dims[1]} & {dims[2]}")))

  interactive_plots <- map(plots, ~ggplotly(.x, tooltip = "text"))

  return(interactive_plots)
}

#' Renders a markdown template file, filtering for a given tissue and outliers
#'
#' @param markdown_path The absolute path to the markdown fle
#' @param tissue Character vector of length 1, indicating the tissue to keep,
#' as present in the metadata sheet for the project
#' @param exclude Character vector of length 1, indicating the Sample IDS to
#' remove seperated by commas, as present in the metadata sheet for the project
#'
#' @return Returns nothing, but creates the corresponding .html files with descriptive names
#' @export
render_markdownfile_single_tisue_updated <- function (markdown_path, tissue, exclude, name=""){

  script_name <- markdown_path %>%
    basename() %>%
    tools::file_path_sans_ext()

  date <- format(Sys.time(), "%Y-%m-%d-%H.%M")

  report_name <- stringr::str_glue("{date}_{script_name}_{tissue}_{name}")

  rmarkdown::render(markdown_path,
                    output_file = report_name,
                    params = list(tissue = tissue,
                                  exclude = exclude,
                                  report_name = report_name),
                    envir = parent.frame())
}
