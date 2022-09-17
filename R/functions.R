# Options -------------------------------------------------------------------------
theme_set(theme_bw())
only_proj_45 <- T

# Constants ---------------------------------------------------------------
factor_cols <- list("sex"=c("1"="male",
                            "2"="female"),
                    "insres_status"=c("1"="res",
                                      "2"="non_res"))

interesting_cols <- c("bmi")


# Functions ---------------------------------------------------------------
## Load data -------------------------------------------------------------------------

#' Load RNA counts from the .rds file
#'
#' @return Count matrix for each sample with NCBI gene names
#' @export
load_rna_counts <- function() {

  g_counts_raw <- readRDS(here("data-raw/salmon.merged.gene_counts.rds"))
  rna_counts <- g_counts_raw %>%
    SummarizedExperiment::assays() %>%
    .$counts

  return(rna_counts)}


#' Construct a corresponding meta data tibble from the salmon count matrix
#'
#' @param rna_counts
#'
#' @return Meta data tibble with Sample.ID, Family, Family-num, Tissue and type.
#' @export
make_salmon_metadata <- function() {
  rna_counts <- load_rna_counts()
  mixed_ids <- c("F52_4_fat",
                 "F52_4_muscle",
                 "F52_4_muscle_dupl",
                 "F53_12_muscle_dupl")

  colnames(rna_counts) %>%
    str_replace("X","") %>%
    str_replace("\\.", "_") %>%
    enframe(name = "Sample.ID", value="Salmon.ID") %>%
    separate(Salmon.ID, into=c("Family", "Family_Num","Tissue", "Type"), remove = F) %>%
    unite(col =sourceId, Family, Family_Num, sep = "-", remove = F) %>%
    replace_na(list(Type="single")) %>%
    mutate(Tissue_fixed=case_when(Salmon.ID %in% mixed_ids & Tissue== "fat" ~ "muscle",
                                  Salmon.ID %in% mixed_ids & Tissue== "muscle" ~ "fat",
                                  TRUE ~ Tissue))
}

#' Load Phnomics platform long phenotypic data file
#'
#' @return A tibble of the long phenotypic data
#' @export
#'
get_pheno_data <- function(){
  if(only_proj_45) {
    df <- read_delim("data-raw/FamilyStudy_PhenoData_13092022.csv") %>%
      filter(projectId==45)}
  else {df <- read_delim("data-raw/FamilyStudy_PhenoData_13092022.csv")}
}
#' Loads and converts the long phenomics metadata to wide format, fixing datatypes
#'
#' @return A wide, deudpped phenomics platform tibble
#' @export
get_pheno_data_wide <- function() {
  pheno_data <- get_pheno_data()

  formatted_wide <- pheno_data %>%
    pivot_wider(names_from = phenoVariable,
                values_from = result,
                id_cols = c(sourceId)) %>%
    mutate(across(3:last_col(), ~str_replace(.x, ",", "."))) %>% # change to english decimal
    mutate(across(3:last_col(),~as.numeric(.x))) # convert to numeric cols %>%


  de_dupped <- formatted_wide %>%
    group_by(sourceId) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm=T)))

  # convert into factors with encoding corresponding to factor_cols
  factored <- de_dupped %>%
    mutate(across(all_of(names(factor_cols)),
                  .fns = function(column) {
                    factor_cols[[cur_column()]][as.character(column)] %>%
                      as.factor()}))

  return(factored)
}

#' Get a phenotypic meta data table for the RNA seq analysis
#'
#' @return Wide, dedupped data for each sample that we have RNAseq data for
#' @export
#'
get_combined_pheno_data <- function(){
  meta_data <- make_salmon_metadata()
  pheno_data <- get_pheno_data_wide()

  combined <- inner_join(meta_data,pheno_data)
  return(combined)
}

#' Load clinical data from DNA methylation paper
#'
#' @return A tibble of the wide metadata
#' @export
#'
get_clinical_data <- function(){
  df <- readxl::read_xls("data-raw/Family_project_data_clinical_characteristics_L.G..xls")
}

# Misc --------------------------------------------------------------------



#' Search for phenotype and get counts in the dataset
#'
#' @param my_string Part of the phenotype to search for
#'
#' @return Returns the phenotype if it exists a long with the number of counts
#' @export

search_for_phenotype <- function(my_string) {

  phenomics_meta_data <- get_combined_pheno_data()

  phenotype_counts <- phenomics_meta_data %>%
    summarise(across(.cols = everything(),
                     .fns = ~.x %>% is.na() %>% `!` %>%  sum())) %>%
    as.vector() %>%
    unlist() %>%
    base::sort(decreasing = T)

  phenotype_counts %>%
    setNames(names(.), .) %>%
    str_detect(my_string) %>%
    phenotype_counts[.]
}


#' Get top 50 most common phenotypes
#'
#' @return The top 50 most common phenotypes along with the counts
#' @export
top_50_most_common_phenotypes <- function() {

  phenomics_meta_data <- get_combined_pheno_data()

  phenotype_counts <- phenomics_meta_data %>%
    summarise(across(.cols = everything(),
                     .fns = ~.x %>% is.na() %>% `!` %>%  sum())) %>%
    as.vector() %>%
    unlist() %>%
    base::sort(decreasing = T)

  phenotype_counts
}



# Deprecated ------------------------------------------------------------------------
