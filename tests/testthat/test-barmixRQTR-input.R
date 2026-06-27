test_that("barmixRQTR checks required inputs before fitting Stan models", {
  valid_data <- list(
    data_count = matrix(
      c(10, 12,
        11, 13,
        14, 16,
        15, 17),
      nrow = 4,
      byrow = TRUE
    ),
    condition_count = c("DMSO", "DMSO", "TreatA", "TreatA"),
    V = c(100, 110, 80, 90),
    condition_v = c("DMSO", "DMSO", "TreatA", "TreatA"),
    cell_line = c("clone1", "clone2")
  )
  
  expect_error(
    barmixRQTR(data = valid_data),
    "`time_d` is required and must be a finite numeric scalar.",
    fixed = TRUE
  )
  
  expect_error(
    barmixRQTR(data = NULL, time_d = 10),
    "`data` must be a list.",
    fixed = TRUE
  )
  
  bad_data <- valid_data
  bad_data$data_count <- NULL
  
  expect_error(
    barmixRQTR(data = bad_data, time_d = 10),
    "Missing required fields in `data`: data_count",
    fixed = TRUE
  )
  
  bad_data <- valid_data
  bad_data$data_count[1, 1] <- NA
  
  expect_error(
    barmixRQTR(data = bad_data, time_d = 10),
    "`data$data_count` contains NA values.",
    fixed = TRUE
  )
  
  bad_data <- valid_data
  bad_data$data_count[1, 1] <- -1
  
  expect_error(
    barmixRQTR(data = bad_data, time_d = 10),
    "`data$data_count` must be non-negative.",
    fixed = TRUE
  )
  
  bad_data <- valid_data
  bad_data$condition_count <- c("DMSO", "TreatA")
  
  expect_error(
    barmixRQTR(data = bad_data, time_d = 10),
    "`data$condition_count` must have the same length as nrow(data$data_count).",
    fixed = TRUE
  )
  
  bad_data <- valid_data
  bad_data$cell_line <- c("clone1")
  
  expect_error(
    barmixRQTR(data = bad_data, time_d = 10),
    "`data$cell_line` must have the same length as ncol(data$data_count).",
    fixed = TRUE
  )
  
  expect_error(
    barmixRQTR(
      data = valid_data,
      time_d = 10,
      control = list(chains = 0)
    ),
    "`control$chains` must be a positive integer.",
    fixed = TRUE
  )
  
  expect_error(
    barmixRQTR(
      data = valid_data,
      time_d = 10,
      dispersion = list(psi_sd = -1)
    ),
    "`dispersion$psi_sd` must be > 0.",
    fixed = TRUE
  )
  
  expect_error(
    barmixRQTR(
      data = valid_data,
      time_d = 10,
      control_group = "missing_control"
    ),
    "`control_group` must be present in `data$condition_count`.",
    fixed = TRUE
  )
})