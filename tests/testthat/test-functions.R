test_that("search_for_phenotype",{
  output <-search_for_phenotype("BMI")
  expect_vector(output, ptype = integer())})


# test_that("get_pheno_data_wide", {get_pheno_data_wide()})

