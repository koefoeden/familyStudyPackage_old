
# Imports -----------------------------------------------------------------
#' @import org.Hs.eg.db
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
#' @import eulerr
#' @import cowplot
NULL


# Options -------------------------------------------------------------------------
library(org.Hs.eg.db)
theme_set(theme_bw())
here::i_am("DESCRIPTION")
only_proj_45 <- T
options(na.action = "na.pass")



# Constants ---------------------------------------------------------------
factor_cols <- list(
  "sex" = c(
    "1" = "male",
    "2" = "female"
  ),
  "insres_status" = c(
    "1" = "res",
    "2" = "non_res"
  )
)

tissues <- c("muscle","fat")
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

  return(rna_counts)
}


#' Construct a corresponding meta data tibble from the salmon count matrix
#' @return Meta data tibble with Sample.ID, Family, Family-num, Tissue and type.
#' @export
make_salmon_metadata <- function() {
  rna_counts <- get_rna_counts()
  mixed_ids <- c(
    "F52_4_fat",
    "F52_4_muscle",
    "F52_4_muscle_dupl",
    "F53_12_muscle_dupl"
  )

  colnames(rna_counts) %>%
    str_replace("X", "") %>%
    str_replace("\\.", "_") %>%
    enframe(name = "Sample.ID", value = "Salmon.ID.edited") %>%
    separate(Salmon.ID.edited, into = c("Family", "Family_Num", "Tissue", "Type"), remove = F) %>%
    unite(col = sourceId, Family, Family_Num, sep = "-", remove = F) %>%
    replace_na(list(Type = "single")) %>%
    mutate(Tissue_fixed = case_when(
      Salmon.ID.edited %in% mixed_ids & Tissue == "fat" ~ "muscle",
      Salmon.ID.edited %in% mixed_ids & Tissue == "muscle" ~ "fat",
      TRUE ~ Tissue
    )) %>%
    mutate(Salmon.ID.raw = colnames(rna_counts))
}

#' Load Phnomics platform long phenotypic data file
#'
#' @return A tibble of the long phenotypic data
#' @export
#'
get_pheno_data <- function() {
  if (only_proj_45) {
    df <- read_delim("data-raw/FamilyStudy_PhenoData_13092022.csv") %>%
      filter(projectId == 45)
  } else {
    df <- read_delim("data-raw/FamilyStudy_PhenoData_13092022.csv")
  }
}

#' Loads and converts the long phenomics metadata to wide format, fixing datatypes
#'
#' @return A wide, deudpped phenomics platform tibble
#' @export
get_pheno_data_wide <- function() {
  pheno_data <- get_pheno_data()

  formatted_wide <- pheno_data %>%
    pivot_wider(
      names_from = phenoVariable,
      values_from = result,
      id_cols = c(sourceId)
    ) %>%
    mutate(across(3:last_col(), ~ str_replace(.x, ",", "."))) %>% # change to english decimal
    mutate(across(3:last_col(), ~ as.numeric(.x))) # convert to numeric cols %>%

  # convert into factors with encoding corresponding to factor_cols

  return(formatted_wide)
}


#' Get a phenotypic meta data table for the RNA seq analysis
#'
#' @return Wide, dedupped data for each sample that we have RNAseq data for
#' @export
#'
get_combined_pheno_data <- function(type = "inner", pheno_data_type = "clinical") {
  if (pheno_data_type == "clinical") {
    pheno_data <- get_clinical_data()
  }

  if (pheno_data_type == "phenomics") {
    pheno_data <- get_pheno_data_wide()
  }

  meta_data <- make_salmon_metadata()

  if (type == "inner") {
    combined <- inner_join(meta_data, pheno_data)
  } else if (type == "left") {
    combined <- left_join(meta_data, pheno_data)
  }
  return(combined)
}

#' Load clinical data from DNA methylation paper. Performs a wide range of data formatting steps
#' and adds some new computed variables such as disease risk based on waist-hip ratio and
#' BIGG_SI measure adjusted for WH
#'
#' @return A tibble of the wide metadata
#' @export
#'
get_clinical_data <- function() {
  mixed_sexes <- c("609-7")
  df <- read_xlsx(here("data-raw/clinical_data.xlsx")) %>%
    mutate(
      SEX = recode_factor(SEX, "1" = "Male", "2" = "Female"), # convert to factors
      smoker = recode_factor(smoker, "2" = "No", "1" = "Yes"),
      Glu_tol_3_class = as.factor(Glu_tol_3_class),
      FAMILY_NR_new = as.factor(FAMILY_NR_new),
      `Far (1) mor (2) syskon (3)barn (4) med T2D` = as.factor(`Far (1) mor (2) syskon (3)barn (4) med T2D`),
      FAMILY_NR_old = as.factor(FAMILY_NR_old),
      Glu_tol_5_class = as.factor(Glu_tol_5_class),
      DATE_LAB = as.Date.character(DATE_LAB, format = "%y%m%d")
    ) %>%
    mutate(WH_risk = case_when(
      WH <= 0.95 & SEX == "Male" ~ "Low",
      WH %>% between(0.96, 1.0) & SEX == "Male" ~ "Moderate",
      WH > 1.0 & SEX == "Male" ~ "High",
      WH <= 0.8 & SEX == "Female" ~ "Low",
      WH %>% between(0.81, 0.85) & SEX == "Female" ~ "Moderate",
      WH > 0.85 & SEX == "Female" ~ "High",
      TRUE ~ "Not classified"
    ))

  data_w_adjusted <- map(
    c("Female", "Male"),
    ~ {
      data <- df %>%
        filter(
          SEX == .x,
          !is.na(BIGG_SI_t_0_60_120)
        )

      model <- glm(
        formula = BIGG_SI_t_0_60_120 ~ WH,
        data = data
      )

      data_w_adjusted <- data %>%
        mutate(BIGG_SI_t_0_60_120_WH_adjusted = model$residuals)

      return(data_w_adjusted)
    }
  )

  combined_data <- bind_rows(
    data_w_adjusted,
    filter(df, is.na(BIGG_SI_t_0_60_120))
  )

  return(combined_data)
}

get_pheno_data_w_RNA <- function(results_dir) {
  combined_pheno_data <- get_combined_pheno_data()

  y_all <- readRDS(file.path(results_dir, "y_all.RDS"))

  cpm_df_all <- cpm(y_all, log = F)

  cpm_transposed <- cbind(nms = names(cpm_df_all), t(cpm_df_all)) %>%
    as_tibble(rownames = "Salmon.ID.raw")

  meta_data_sheet_w_RNA_seq <-
    left_join(
      x = cpm_transposed,
      y = combined_pheno_data,
      by = "Salmon.ID.raw",
    )

  return(meta_data_sheet_w_RNA_seq)
}

#' Get PCA data from Plink PCA output by loading file and formatting it
#'
#' @param eigenvec_data The path to the PLINK PCA output
#'
#' @return Formatted PCA data with sourceId as identifying column
#' @export
get_PCA_data <- function(eigenvec_data = "data-raw/fam_study_only_genotyped.eigenvec") {
  PCA_data_plink <- read_tsv(file = eigenvec_data)

  PCA_data <- PCA_data_plink %>%
    mutate(originalID = `#IID`) %>%
    relocate(originalID) %>%
    filter(str_detect(`#IID`, pattern = "41x", negate = T)) %>% # remove non 45 samples
    mutate(`#IID` = str_remove_all(`#IID`, "45x")) %>% # remove 45 designation
    separate(`#IID`, into = c("fam", "ind"), sep = c("-")) %>% # Seperate IID into fam and individuals
    separate(fam, into = c("#IID", "fam"), sep = c("_")) %>% # Separate IID into IID and family
    # extract(ind, into=c("ind", "replicate"), regex = "([0-9]{1,2}):?([A-Z])?") %>%  # add into replicates
    # mutate(replicate=na_if(replicate, "")) %>%
    mutate(fam = coalesce(fam, `#IID`)) %>% # if fam is empty, take value fro IID
    filter(!is.na(ind)) %>%
    mutate(sourceId = paste(fam, ind, sep = "-")) %>% # remove rows with missing info
    relocate(sourceId)

  return(PCA_data)
}

#' Get pheno data w genetic info
#'
#' @param results_dir Path to a RNA seq analysis data folder
#' @return Returns a tibble of all samples that have RNA data with their
#' clinical characteristics, RNA and genetic data
#' @export
get_pheno_data_w_genetics <- function(results_dir) {
  pheno_data_w_RNA <- get_pheno_data_w_RNA(results_dir = results_dir)
  PCA_data <- get_PCA_data()

  pheno_data_w_genetics <- left_join(x = pheno_data_w_RNA, y = PCA_data)
  return(pheno_data_w_genetics)
}

#' Writes a gene info data frame to RDS with the ENSEMBL ID, HGNC symbol and gene description
#'
#' @param genes_name_type "NCBI" or "ENSEMBL"
#' @param organism "human", "rat" or "mouse"
#'
#' @return Nothing.
#' @export
write_gene_info_df <- function(genes_name_type="NCBI", organism="human") {

  # get the salmon RNA counts object to figure out which genes to get info for
  rna_counts <- get_rna_counts()

  # get the mart containing the data
  mart <-  readRDS(file = if (organism=="mouse") here("shared_data/mouse_mart.RDS")
                   else if (organism=="human") here("shared_data/human_mart.RDS")
                   else if (organism=="rat") here("shared_data/rat_mart.RDS"))

  # Get the relevant info for each rowname in counts
  gene_info_df<- getBM(filters= case_when(genes_name_type=="NCBI" ~"entrezgene_accession",
                                          genes_name_type=="ENSEMBL" ~ "ensembl_gene_id"),
                       attributes= c("entrezgene_accession", "ensembl_gene_id","description"),
                       values=rownames(rna_counts),
                       mart=mart) %>%
    # remove 1:2 mappings
    distinct(across(case_when(genes_name_type=="NCBI" ~"entrezgene_accession",
                              genes_name_type=="ENSEMBL" ~ "ensembl_gene_id")),
             .keep_all = TRUE) %>%
    # rename for better understanding
    dplyr::rename("Name"=entrezgene_accession,
                  "ENSEMBL_ID"=ensembl_gene_id) %>%
    # format description
    mutate(description=str_replace(string = description,
                                   pattern = " \\[Source:.*\\]",
                                   replacement = ""))


  saveRDS(object = gene_info_df,
          file = str_glue("shared_data/gene_info_df_{genes_name_type}_centric.RDS"))

}

#' Creates a fitted DGE list from tissue vector and model formula stirng
#'
#' @param tissue the tissue to keept
#' @param formula the formula for the model as a string
#'
#' @return the fitted DGE list
#' @export
get_whole_DGE_list_by_tissue <- function(tissue="muscle") {

  # get metadata table for each sample
  combined_pheno_data <- get_combined_pheno_data()
  combined_pheno_data_single_tissue <- combined_pheno_data %>%
    filter(Tissue_fixed == tissue) %>%
    filter(sourceId!="641-9")

  ## DGE list creation for all samples
  rna_counts_unordered <- get_rna_counts()
  rna_counts <- rna_counts_unordered[,as.character(combined_pheno_data_single_tissue$Salmon.ID.raw)]
  y_all <- DGEList(rna_counts,
                   samples=combined_pheno_data_single_tissue)

  # add CPM counts to samples slot
  cpm_df_all <- cpm(y_all, log = F)

  cpm_transposed <- cbind(nms = names(cpm_df_all), t(cpm_df_all)) %>%
    as_tibble(rownames = "Salmon.ID.raw")

  samples_info_w_cpm <-
    left_join(
      x = y_all$samples,
      y = cpm_transposed,
      by = "Salmon.ID.raw",
    )

  y_all$samples <- samples_info_w_cpm

  return(y_all)
}

get_subsetted_DGE_list_by_tissue_formula <- function(tissue="muscle", formula="~BIGG_SI_t_0_60_120 + SEX + AGE") {

  unfitted_DGE_list_by_tissue <- get_whole_DGE_list_by_tissue(tissue = tissue)

  design_all <-  model.matrix(object = as.formula(formula),
                              data= unfitted_DGE_list_by_tissue$samples,
                              na.action="na.pass")

  colnames(design_all) <- make.names(colnames(design_all))

  idx_expressed_all <- filterByExpr(unfitted_DGE_list_by_tissue, design = design_all)


  subsetted_DGE_list <- unfitted_DGE_list_by_tissue[idx_expressed_all,] %>%
    calcNormFactors() %>%
    estimateDisp(design=design_all)

  return(subsetted_DGE_list)
}

get_coefficients_from_DGE_list <- function(DGE_list) {
  # get all coefficients in linear model except intercept
  coefficients <- DGE_list$design %>%
    colnames() %>%
    `[`(-1) %>%
    set_names()

  return(coefficients)
}
#' Performs the test for each coefficient in the linear model
#'
#' @param fit_all the fitted DGE list
#' @param gene_names_type "NCBI" or "ENSEMBL"
#'
#' @return
#' @export
get_DGE_data <- function(y_all, gene_names_type="NCBI") {

  # fit data
  fit_all <- glmQLFit(y_all, robust = TRUE)

  coefficients <- get_coefficients_from_DGE_list(fit_all)

  DE_genes_test_all <- coefficients %>%
    purrr::map(~edgeR_tester(coef = .x,
                             efit = fit_all,
                             id_name = "Name"))

  gene_info_df <- readRDS(str_glue("shared_data/gene_info_df_{gene_names_type}_centric.RDS"))

  DE_genes_test_all_w_info <- DE_genes_test_all %>%
    map(~dplyr::full_join(., y=gene_info_df))

  return(DE_genes_test_all_w_info)
}

get_GO_data <- function(fit_all, gene_names_type="NCBI", organism="human") {


  GO_keytype <- case_when(gene_names_type=="ENSEMBL" ~ "ENSEMBL",
                          gene_names_type=="NCBI" ~ "SYMBOL")

  # GO_terms_all <- readRDS(str_glue("shared_data/GO_terms_BP_{organism}_{GO_keytype}.RDS"))
  # not equivalent to below??

  GO_terms_all <- get_enrichment_terms(org_db =org.Hs.eg.db,
                                       gene_ids = rownames(fit_all),
                                       gene_id_key_type = "SYMBOL",
                                       min_genes = 5, max_genes = 500)$BP


  coefficients <- get_coefficients_from_DGE_list(fit_all)

  ontology_tests_all <- run_ontology_tests(GO_terms_all,
                                           fit_all,
                                           coefs = coefficients,
                                           fun=limma::camera)

  return(ontology_tests_all)
}

run_full_analysis <- function(tissues=c("muscle","fat"),
                              formulas ="~BIGG_SI_t_0_60_120 + SEX + AGE") {
  output <- list()
  for (tissue in tissues) {
    output[[tissue]] <- list()
    whole_DGE_list <- get_whole_DGE_list_by_tissue(tissue = tissue)

    output[[tissue]][["DATA"]] <- whole_DGE_list

    for (formula in formulas) {

      output[[tissue]][[formula]] <- list()
      subsetted_DGE_list <- get_subsetted_DGE_list_by_tissue_formula(tissue = tissue,
                                                                     formula = formula)

      DGE_all <- get_DGE_data(subsetted_DGE_list)
      output[[tissue]][[formula]][["DGE"]] <- DGE_all

      GO_all <- get_GO_data(subsetted_DGE_list)
      output[[tissue]][[formula]][["GO"]] <- GO_all


    }
  }
  time <- format(Sys.time(), "%d_%H_%M")
  tissues_formatted <- paste(tissues, collapse = "_")
  formulas_formatted <- str_replace(formulas, "~", "") %>%
    str_replace_all(" \\+ ", "_")

  save_string <- str_glue("results_{tissues_formatted}_{formulas_formatted}_{time}.RDS")
  saveRDS(output, save_string)
  cat("Saved results to:", save_string)

  return(output)
}
# Misc --------------------------------------------------------------------

#' Get interesting genes
#'
#' @return Vector of gene names
#' @export
get_interesting_genes <- function() {
  gene_vector <- read_xlsx(here("data-raw/interesting_genes.xlsx")) %>%
    pull(1) %>%
    na.omit()
  return(gene_vector)
}


#' Function to get genes from a table https://www-nature-com.ep.fjernadgang.kb.dk/articles/s41467-019-13869-w#Sec1
#' that includes exercise induced genes. Filtering for the long term
#'
#' @return The chararcter vector of genes
get_long_term_exercise_induced_genes <- function(){

  genes <- read_xlsx("data-raw/exercise_genes.xlsx") %>%
    mutate(across(2:last_col(), #format to remove weird spaces and minus sign and convert to numeric
                  ~.x %>% str_replace(pattern = "−", replacement = "-") %>%
                    gsub(pattern = "\\s", replacement = "") %>%
                    as.numeric())) %>%
    filter(`Training aerobic, FDR`<0.05 | `Training resistance, FDR`<0.05) %>%
    pull(NAME)

  return(genes)
}

get_SLC2A4_TFs <- function() {
  gene_vector <- read_xlsx(here("data-raw/SLC2A4_TFs.xlsx")) %>%
    pull(1) %>%
    unique()
  return(gene_vector)
}
#' Get dataframe of interesting genes
#'
#' @return Vector of gene names
#' @export
get_interesting_genes_df <- function() {
  gene_df <- read_xlsx(here("data-raw/interesting_genes.xlsx"))
  return(gene_df)
}


#' Search for phenotype and get counts in the dataset
#'
#' @param my_string Part of the phenotype to search for
#'
#' @return Returns the phenotype if it exists a long with the number of counts
#' @export

search_for_phenotype <- function(my_string) {
  phenomics_meta_data <- get_combined_pheno_data()

  phenotype_counts <- phenomics_meta_data %>%
    summarise(across(
      .cols = everything(),
      .fns = ~ .x %>%
        is.na() %>%
        `!`() %>%
        sum()
    )) %>%
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
    summarise(across(
      .cols = everything(),
      .fns = ~ .x %>%
        is.na() %>%
        `!`() %>%
        sum()
    )) %>%
    as.vector() %>%
    unlist() %>%
    base::sort(decreasing = T)

  phenotype_counts
}

# Plots -------------------------------------------------------------------
plot_PCA <- function(n_PCs = 10) {
  PCA_data <- get_PCA_data()

  plots_for_each_PC <- map(
    seq(1, n_PCs, 2),
    ~ {
      my_plot <- PCA_data %>%
        ggplot(aes_string(
          x = str_glue("PC{.x}"),
          y = str_glue("PC{.x+1}"),
          color = "fam",
          label = "ind"
        )) +
        geom_point() +
        theme_bw() +
        theme(legend.position = "bottom")

      interactive_plot <- ggplotly(my_plot)
      return(interactive_plot)
    }
  )
  return(plots_for_each_PC)
}

plot_oral_glucose_data <- function() {
  combined_pheno_data <- get_combined_pheno_data()

  plots <- map(c("glu", "ins", "cpep"),
               .f = ~ combined_pheno_data %>%
                 group_by(Glu_tol_3_class) %>%
                 summarise(across(starts_with("OGTT"), ~ median(.x, na.rm = T))) %>%
                 pivot_longer(
                   cols = matches(str_glue("OGTT_{.x}_t-*[0-9]+$")),
                   names_prefix = str_glue("OGTT_{.x}_t"),
                   names_to = "Timepoint",
                   values_to = str_glue("{.x}") %>% as.character(),
                   names_transform = as.numeric,
                   values_transform = as.numeric
                 ) %>%
                 ggplot(aes_string(
                   x = "Timepoint",
                   y = str_glue("{.x}"),
                   color = "Glu_tol_3_class"
                 )) +
                 geom_line(alpha = 0.5)
  )
  return(plots)
}

# Deprecated ------------------------------------------------------------------------
