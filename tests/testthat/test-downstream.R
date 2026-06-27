make_mock_barmixR_fit <- function() {
  structure(
    list(
      names_cell = c("clone1", "clone2"),
      condition_count = c("DMSO", "DMSO", "TreatA", "TreatA"),
      condition_v = c("DMSO", "DMSO", "TreatA", "TreatA"),
      control_group = "DMSO",
      n_sam = 20L,
      L = 2L,
      K = 2L,
      time_d = 10,
      in_vivo = TRUE
    ),
    class = "barmixR_fit"
  )
}

test_that("fractionRatio validates inputs", {
  model <- make_mock_barmixR_fit()
  
  expect_error(
    fractionRatio(model = list(), sampled_fraction = list()),
    "`model` must be an object of class 'barmixR_fit' returned by `barmixRQTR()`.",
    fixed = TRUE
  )
  
  expect_error(
    fractionRatio(model = model, sampled_fraction = NULL),
    "`sampled_fraction` must be a list of posterior samples from `ppcBarcodes()`.",
    fixed = TRUE
  )
})

test_that("fractionRatio returns expected output structure", {
  model <- make_mock_barmixR_fit()
  
  sampled_fraction <- list(
    DMSO = matrix(runif(40, 0.1, 1), nrow = 20, ncol = 2),
    TreatA = matrix(runif(40, 0.1, 1), nrow = 20, ncol = 2)
  )
  
  result <- fractionRatio(
    model = model,
    sampled_fraction = sampled_fraction
  )
  
  expect_s3_class(result, "barmixR_result")
  expect_named(result, c("plot_ratio_fraction", "li_sam_ratio_relative", "type"))
  expect_s3_class(result$plot_ratio_fraction, "gg")
  expect_type(result$li_sam_ratio_relative, "list")
  expect_equal(result$type, "fraction_ratio")
})

test_that("populationRatio validates inputs", {
  model <- make_mock_barmixR_fit()
  
  expect_error(
    populationRatio(model = list(), sampled_elements = list()),
    "`model` must be an object of class 'barmixR_fit' returned by `barmixRQTR()`.",
    fixed = TRUE
  )
  
  expect_error(
    populationRatio(model = model, sampled_elements = NULL),
    "`sampled_elements` must be a list of posterior samples from `ppcPopulation()`.",
    fixed = TRUE
  )
})

test_that("populationRatio returns expected output structure", {
  model <- make_mock_barmixR_fit()
  
  sampled_population <- list(
    DMSO = runif(20, 90, 120),
    TreatA = runif(20, 70, 100)
  )
  
  result <- populationRatio(
    model = model,
    sampled_elements = sampled_population
  )
  
  expect_s3_class(result, "barmixR_result")
  expect_named(result, c("plot_ratio_population", "merged_data", "li_sam_ratio_V", "type"))
  expect_s3_class(result$plot_ratio_population, "gg")
  expect_s3_class(result$merged_data, "data.frame")
  expect_type(result$li_sam_ratio_V, "list")
  expect_equal(result$type, "population_ratio")
})

test_that("QTRresistance validates inputs", {
  model <- make_mock_barmixR_fit()
  
  expect_error(
    QTRresistance(model = list(), li_sam_ratio_relative = list(), li_sam_ratio_V = list()),
    "`model` must be an object of class 'barmixR_fit' returned by `barmixRQTR()`.",
    fixed = TRUE
  )
  
  expect_error(
    QTRresistance(model = model, li_sam_ratio_relative = NULL, li_sam_ratio_V = list()),
    "`li_sam_ratio_relative` must be a list of posterior fraction ratios from fractionRatio().",
    fixed = TRUE
  )
  
  expect_error(
    QTRresistance(model = model, li_sam_ratio_relative = list(), li_sam_ratio_V = NULL),
    "`li_sam_ratio_V` must be a list of posterior population ratios from populationRatio().",
    fixed = TRUE
  )
})

test_that("QTRresistance returns expected output structure", {
  model <- make_mock_barmixR_fit()
  
  li_sam_ratio_relative <- list(
    TreatA = matrix(runif(40, 0.5, 1.5), nrow = 20, ncol = 2)
  )
  
  li_sam_ratio_V <- list(
    TreatA = runif(20, 0.5, 1.5)
  )
  
  result <- QTRresistance(
    model = model,
    li_sam_ratio_relative = li_sam_ratio_relative,
    li_sam_ratio_V = li_sam_ratio_V
  )
  
  expect_s3_class(result, "barmixR_result")
  expect_named(result, c("treatment_resistance", "summary_table", "type"))
  expect_s3_class(result$treatment_resistance, "gg")
  expect_s3_class(result$summary_table, "data.frame")
  expect_true(all(c(
    "cell", "treat", "mean", "median", "sd",
    "q25", "q75", "lower_whisker", "upper_whisker", "CDF"
  ) %in% colnames(result$summary_table)))
  expect_equal(result$type, "QTR")
})

test_that("QTRDecision returns expected output structure", {
  model <- make_mock_barmixR_fit()
  
  summary_table <- data.frame(
    cell = rep(c("clone1", "clone2"), each = 2),
    treat = rep(c("TreatA", "TreatB"), times = 2),
    median = c(0.1, -0.2, 0.3, -0.1),
    q25 = c(0.0, -0.3, 0.2, -0.2),
    q75 = c(0.2, -0.1, 0.4, 0.0),
    lower_whisker = c(-0.1, -0.4, 0.1, -0.3),
    upper_whisker = c(0.3, 0.0, 0.5, 0.1)
  )
  
  result <- QTRDecision(
    model = model,
    summary_table = summary_table,
    ncol = 1
  )
  
  expect_s3_class(result, "barmixR_result")
  expect_named(result, c("rank_plot", "rank_summary", "type"))
  expect_s3_class(result$rank_summary, "data.frame")
  expect_true("rank" %in% colnames(result$rank_summary))
  expect_equal(result$type, "QTR_decision")
})

test_that("QTRheatmap returns expected output structure", {
  model <- make_mock_barmixR_fit()
  
  summary_table <- data.frame(
    cell = rep(c("clone1", "clone2"), each = 2),
    treat = rep(c("DMSO", "TreatA"), times = 2),
    mean = c(0.1, -0.2, 0.3, -0.1),
    CDF = c(0.8, 0.2, 0.7, 0.3)
  )
  
  result <- QTRheatmap(
    model = model,
    summary_table = summary_table
  )
  
  expect_s3_class(result, "barmixR_result")
  expect_named(result, c("heatmap_within", "heatmap_treatment", "type"))
  expect_s3_class(result$heatmap_within, "gg")
  expect_s3_class(result$heatmap_treatment, "gg")
  expect_equal(result$type, "QTR_heatmap")
})