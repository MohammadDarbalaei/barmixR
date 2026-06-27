test_that("ppcBarcodes validates model class", {
  expect_error(
    ppcBarcodes(model = list()),
    "`model` must be an object of class 'barmixR_fit' returned by `barmixRQTR()`.",
    fixed = TRUE
  )
})

test_that("ppcPopulation validates model class and required fields", {
  expect_error(
    ppcPopulation(model = list()),
    "`model` must be an object of class 'barmixR_fit' returned by `barmixRQTR()`.",
    fixed = TRUE
  )
  
  incomplete_model <- structure(
    list(condition_v = c("DMSO", "TreatA")),
    class = "barmixR_fit"
  )
  
  expect_error(
    ppcPopulation(model = incomplete_model),
    "Missing required fields in `model`: fit_V, V, n_sam",
    fixed = TRUE
  )
})