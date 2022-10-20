
# Imports -----------------------------------------------------------------
#' @import ggplot2
#' @import dplyr
#' @import edgeR
#' @import purrr
#' @import readr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import forcats
#' @import DT
#' @import SummarizedExperiment
#' @import cbmr
#' @import gtools
#' @import readxl
#' @import here
#' @import plotly
#' @import ggVennDiagram
NULL


# Options -------------------------------------------------------------------------
theme_set(theme_bw())
here::i_am("DESCRIPTION")
only_proj_45 <- T
options(na.action='na.pass')


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
get_rna_counts <- function() {

  g_counts_raw <- readRDS(here("data-raw/salmon.merged.gene_counts.rds"))
  rna_counts <- g_counts_raw %>%
    SummarizedExperiment::assays() %>%
    .$counts

  return(rna_counts)}


#' Construct a corresponding meta data tibble from the salmon count matrix
#' @return Meta data tibble with Sample.ID, Family, Family-num, Tissue and type.
#' @export
make_salmon_metadata <- function() {
  rna_counts <- get_rna_counts()
  mixed_ids <- c("F52_4_fat",
                 "F52_4_muscle",
                 "F52_4_muscle_dupl",
                 "F53_12_muscle_dupl")


  colnames(rna_counts) %>%
    str_replace("X","") %>%
    str_replace("\\.", "_") %>%
    enframe(name = "Sample.ID", value="Salmon.ID.edited") %>%
    separate(Salmon.ID.edited, into=c("Family", "Family_Num","Tissue", "Type"), remove = F) %>%
    unite(col = sourceId, Family, Family_Num, sep = "-", remove = F) %>%
    replace_na(list(Type="single")) %>%
    mutate(Tissue_fixed=case_when(Salmon.ID.edited %in% mixed_ids & Tissue== "fat" ~ "muscle",
                                  Salmon.ID.edited %in% mixed_ids & Tissue== "muscle" ~ "fat",
                                  TRUE ~ Tissue)) %>%
    mutate(Salmon.ID.raw=colnames(rna_counts))
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

  # convert into factors with encoding corresponding to factor_cols

  return(formatted_wide)
}


#' Get a phenotypic meta data table for the RNA seq analysis
#'
#' @return Wide, dedupped data for each sample that we have RNAseq data for
#' @export
#'
get_combined_pheno_data <- function(type="inner", pheno_data_type="clinical"){
  if (pheno_data_type=="clinical") {
    pheno_data<- get_clinical_data()
  }

  if (pheno_data_type=="phenomics") {
    pheno_data<-get_pheno_data_wide()
  }

  meta_data <- make_salmon_metadata()

  if (type=="inner") {
    combined <- inner_join(meta_data,pheno_data)
  } else if (type=="left") {
    combined <- left_join(meta_data,pheno_data)
  }
  return(combined)
}

#' Load clinical data from DNA methylation paper
#'
#' @return A tibble of the wide metadata
#' @export
#'
get_clinical_data <- function(){
  mixed_sexes <- c("609-7")
  df <- read_xlsx(here("data-raw/clinical_data.xlsx")) %>%
    mutate(SEX=recode_factor(SEX, "1"="Male", "2"="Female"), # convert to factors
           smoker=recode_factor(smoker, "2"="No", "1"="Yes"),
           Glu_tol_3_class=as.factor(Glu_tol_3_class),
           FAMILY_NR_new=as.factor(FAMILY_NR_new),
           `Far (1) mor (2) syskon (3)barn (4) med T2D`=as.factor(`Far (1) mor (2) syskon (3)barn (4) med T2D`),
           FAMILY_NR_old=as.factor(FAMILY_NR_old),
           Glu_tol_5_class=as.factor(Glu_tol_5_class),
           DATE_LAB=as.Date.character(DATE_LAB, format = "%y%m%d")) %>%
    mutate(WH_risk=case_when(WH <= 0.95 & SEX=="Male" ~ "Low",
                             WH %>% between(0.96,1.0) & SEX=="Male" ~ "Moderate",
                             WH > 1.0 & SEX == "Male" ~ "High",
                             WH <= 0.8 & SEX == "Female" ~ "Low",
                             WH %>% between(0.81,0.85) & SEX=="Female" ~ "Moderate",
                             WH > 0.85 & SEX == "Female" ~ "High",
                             TRUE ~ "Not classified"))
  return(df)
}
get_pheno_data_w_RNA <- function(results_dir) {
  combined_pheno_data <- get_combined_pheno_data()

  y_all <- readRDS(file.path(results_dir, "y_all.RDS"))

  cpm_df_all <- cpm(y_all, log=F)

  cpm_transposed <- cbind(nms = names(cpm_df_all), t(cpm_df_all)) %>%
    as_tibble(rownames = "Salmon.ID.raw")

  meta_data_sheet_w_RNA_seq <-
    left_join(x=cpm_transposed,
              y=combined_pheno_data,
              by="Salmon.ID.raw", )

  return(meta_data_sheet_w_RNA_seq)
}

#' Get PCA data from Plink PCA output by loading file and formatting it
#'
#' @param eigenvec_data The path to the PLINK PCA output
#'
#' @return Formatted PCA data with sourceId as identifying column
#' @export
#'
#' @examples
get_PCA_data <- function(eigenvec_data="data-raw/fam_study_only_genotyped.eigenvec") {

  PCA_data_plink <- read_tsv(file = eigenvec_data)


  PCA_data <- PCA_data_plink %>%
    mutate(originalID=`#IID`) %>%
    relocate(originalID) %>%
    filter(str_detect(`#IID`, pattern = "41x", negate = T)) %>% # remove non 45 samples
    mutate(`#IID`=str_remove_all(`#IID`, "45x")) %>%  # remove 45 designation
    separate(`#IID`, into=c("fam","ind"), sep=c("-")) %>% # Seperate IID into fam and individuals
    separate(fam, into=c("#IID","fam"), sep=c("_")) %>% # Separate IID into IID and family
    mutate(fam=coalesce(fam,`#IID`)) %>% # if fam is empty, take value fro IID
    filter(!is.na(ind)) %>%
    mutate(sourceId=paste(fam, ind, sep="-")) %>% # remove rows with missing info
    relocate(sourceId)

  return(PCA_data)
}

#' Get pheno data w genetic info
#'
#' @param results_dir Path to a RNA seq analysis data folder
#' @return Returns a tibble of all samples that have RNA data with their
#' clinical characteristics, RNA and genetic data
#'
#' @export
#'
#' @examples
get_pheno_data_w_genetics <- function(results_dir) {

  pheno_data_w_RNA <- get_pheno_data_w_RNA(results_dir =results_dir )
  PCA_data <- get_PCA_data()

  pheno_data_w_genetics <- left_join(x = pheno_data_w_RNA, y=PCA_data)
  return(pheno_data_w_genetics)
}
# Misc --------------------------------------------------------------------

#' Get interesting genes
#'
#' @return Vector of gene names
#' @export
get_interesting_genes <- function() {read_xlsx(here("interesting_genes.xlsx")) %>% pull(Name)}

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


# Plots -------------------------------------------------------------------
plot_PCA <- function(n_PCs=10) {
  PCA_data <- get_PCA_data()

  plots_for_each_PC <- map(seq(1,n_PCs,2),
                           ~{my_plot <- PCA_data %>%
                             ggplot(aes_string(x=str_glue("PC{.x}"),
                                               y=str_glue("PC{.x+1}"),
                                               color="fam",
                                               label="ind")) +
                             geom_point() +
                             theme_bw() +
                             theme(legend.position = "bottom")

                           interactive_plot <- ggplotly(my_plot)
                           return(interactive_plot)
                           })
  return(plots_for_each_PC)
}

plot_oral_glucose_data <- function() {

  combined_pheno_data <- get_combined_pheno_data()

  plots <- map(c("glu","ins","cpep"),
               .f = ~combined_pheno_data %>%
                 group_by(Glu_tol_3_class) %>%
                 summarise(across(starts_with("OGTT"), ~median(.x, na.rm = T))) %>%
                 pivot_longer(cols=matches(str_glue("OGTT_{.x}_t-*[0-9]+$")),
                              names_prefix = str_glue("OGTT_{.x}_t"),
                              names_to = "Timepoint",
                              values_to = str_glue("{.x}") %>% as.character(),
                              names_transform = as.numeric,
                              values_transform = as.numeric) %>%
                 ggplot(aes_string(x="Timepoint",
                                   y=str_glue("{.x}"),
                                   color="Glu_tol_3_class")) +
                 geom_line(alpha=0.5))
  return(plots)
}

# Deprecated ------------------------------------------------------------------------
