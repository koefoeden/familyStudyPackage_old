% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/single_cell_functions.R
\name{create_multimodal_object}
\alias{create_multimodal_object}
\title{Create and save multimodal Seurat object}
\usage{
create_multimodal_object(
  filtered_h5_feature_bc_matrix_path = "data-raw/filtered_feature_bc_matrix.h5",
  frag_path = "data-raw/atac_fragments.tsv.gz",
  qs_save_path = "data/seur.qs"
)
}
\arguments{
\item{fragpath}{}
}
\value{
Object (and saves it as well)
}
\description{
Create and save multimodal Seurat object
}
