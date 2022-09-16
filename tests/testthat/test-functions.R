setwd(here())

test_that("search_for_phenotype",{
  expect_equal(search_for_phenotype("lep"),
               c("leptin"=231))
})


?test_that

