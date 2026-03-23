#' Posterior predictive checks for population size
#'
#' Performs posterior predictive checks (PPC) for population size data
#' modeled in the barmixR framework. Population size corresponds to
#' tumor volume in \emph{in vivo} experiments or cellular confluency
#' in \emph{in vitro} assays.
#'
#' Posterior predictive samples are generated from the population size
#' component of the barmixR model fitted by \code{barmixRQTR()}. Observed
#' values are compared with predicted values to assess model fit and
#' uncertainty. Predictive distributions are visualized using violin
#' plots with observed measurements overlaid.
#'
#' @param model An object of class \code{\link{barmixR_fit}} returned by \code{\link{barmixRQTR}}
#'   containing:
#'   \itemize{
#'     \item \code{condition_v}: Vector of treatment labels.
#'     \item \code{fit_V}: Fitted population size model
#'       (tumor volume in \emph{in vivo} experiments or
#'       confluency in \emph{in vitro} assays).
#'     \item \code{V}: Observed population size measurements.
#'     \item \code{n_sam}: Number of posterior samples retained.
#'     \item \code{VT0}: Optional baseline population size values used
#'       for normalization.
#'   }
#'
#' @param q_up Upper quantile threshold (unused; kept for API
#'   compatibility).
#'
#' @param q_lo Lower quantile threshold (unused; kept for API
#'   compatibility).
#'
#' @return An object of class \code{barmixR_diagnostic} containing:
#'   \itemize{
#'     \item \code{ppc_plot}: Posterior predictive check plot.
#'
#'     \item \code{ppc_plot_normalized}: Posterior predictive check plot
#'       for normalized population size values (if \code{VT0} is provided).
#'
#'     \item \code{ppc_data}: Data used for visualization.
#'
#'     \item \code{sampled_population}: Posterior predictive samples of
#'       population size.
#'       
#'     \item \code{type}: Result type ("ppc_population"). 
#'   }
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_point labs theme
#' @importFrom ggplot2 element_text element_rect margin scale_y_log10
#'
#' @examples
#' set.seed(1)
#'
#' ## ---------------------------
#' ## Simulate small example data
#' ## ---------------------------
#'
#' # barcode count matrix (4 samples x 5 clones)
#' data_count <- matrix(
#'   c(
#'     20, 22, 24, 23, 21,   # DMSO_rep1
#'     21, 23, 25, 22, 20,   # DMSO_rep2
#'     26, 25, 22, 20, 18,   # TreatA_rep1
#'     25, 24, 23, 21, 19    # TreatA_rep2
#'   ),
#'   nrow = 4,
#'   byrow = TRUE
#' )
#'
#' colnames(data_count) <- paste0("clone", 1:5)
#' rownames(data_count) <- c("DMSO_rep1", "DMSO_rep2",
#'                          "TreatA_rep1", "TreatA_rep2")
#'
#' # treatment labels
#' condition_count <- factor(c("DMSO", "DMSO", "TreatA", "TreatA"))
#'
#' # population size (e.g. tumor volume)
#' V <- c(300, 280, 100, 120)
#'
#' # assemble input
#' data <- list(
#'   data_count = data_count,
#'   condition_count = condition_count,
#'   V = V,
#'   condition_v = condition_count,
#'   cell_line = factor(colnames(data_count))
#' )
#'
#' ## Fit model (fast settings for example only)
#' ## Note: iteration numbers are intentionally small to ensure the example runs quickly.
#' ## For real analysis, increase iterations (e.g., 10000) and use multiple chains.
#'
#' fit <- barmixRQTR(
#'   data = data,
#'   time_d = 10,
#'   control = list(
#'     chains = 1,
#'     iter_count = 50,
#'     iter_V = 50,
#'     cores = 1
#'   ),
#'   seed = 123
#' )
#' 
#' ## ---------------------------
#' ## Posterior Predictive Check for Population Size
#' ## ---------------------------
#'
#' ppc_pop <- ppcPopulation(fit)
#'
#' ## visualize PPC
#' ppc_pop$ppc_plot
#' 
#' @export
ppcPopulation <- function(
    model,
    q_up = 0.975,
    q_lo = 0.025
) {

  ## ---------- helpers ----------
  .stop_if <- function(cond, msg) {
    if (isTRUE(cond)) stop(msg, call. = FALSE)
  }

  ## ---------- validate inputs ----------
  .stop_if(is.null(model) || !inherits(model, "barmixR_fit"),
           "`model` must be an object of class 'barmixR_fit' returned by `barmixRQTR()`.")

  required_fields <- c("condition_v", "fit_V", "V", "n_sam")
  missing_fields <- setdiff(required_fields, names(model))

  .stop_if(length(missing_fields) > 0,
           paste0("Missing required fields in `model`: ",
                  paste(missing_fields, collapse = ", ")))

  .stop_if(!requireNamespace("rstan", quietly = TRUE),
           "Package `rstan` is required but not installed.")

  ## ---------- extract model elements ----------
  condition_v <- as.character(model$condition_v)
  fit_V <- model$fit_V
  V <- as.numeric(model$V)
  n_sam <- as.integer(model$n_sam)
  VT0 <- model$VT0

  conditions <- unique(condition_v)

  ## ---------- extract posterior predictions ----------
  yrep_V <- rstan::extract(fit_V, "ypred")$ypred

  if (is.null(VT0)) {
    uncertainty_V <- yrep_V
  } else {
    uncertainty_V <- yrep_V /
      matrix(rep(VT0, each = nrow(yrep_V)), nrow = nrow(yrep_V))
  }

  ## ---------- group indices ----------
  rep_per_treat <- table(factor(condition_v, levels = conditions))

  start <- 1
  end <- cumsum(rep_per_treat)

  seq_list <- vector("list", length(rep_per_treat))

  for (i in seq_along(rep_per_treat)) {
    seq_list[[i]] <- start:end[i]
    start <- end[i] + 1
  }

  ## ---------- sample predictive uncertainty ----------
  unc_V_treat <- vector("list", length(seq_list))
  unc_V_treat_ppc <- vector("list", length(seq_list))

  for (i in seq_along(seq_list)) {

    unc_V_treat[[i]] <-
      as.vector(uncertainty_V[, seq_list[[i]]])

    unc_V_treat_ppc[[i]] <-
      as.vector(yrep_V[, seq_list[[i]]])
  }

  sampled_elements <-
    lapply(unc_V_treat, function(x) sample(x, n_sam))
  names(sampled_elements) <- conditions
  sampled_elements_ppc <-
    lapply(unc_V_treat_ppc, function(x) sample(x, n_sam))

  ## ---------- prepare plotting data ----------
  gr_lev <- rep(
    factor(conditions, levels = conditions),
    each = n_sam
  )

  approx_data <- data.frame(
    predicted_population = unlist(sampled_elements_ppc),
    treatment = gr_lev,
    stringsAsFactors = FALSE
  )

  approx_data_Nor_V <- data.frame(
    predicted_population = unlist(sampled_elements),
    treatment = gr_lev,
    stringsAsFactors = FALSE
  )

  data_v <- data.frame(
    V = V,
    treat = factor(condition_v, levels = conditions)
  )

  if (!is.null(VT0)) {
    data_Nor_V <- data.frame(
      V = V / VT0,
      treat = factor(condition_v, levels = conditions)
    )
  }

  ## ---------- determine labels ----------
  y_label <- if (all(V >= 0 & V <= 1)) {
    "Confluency & uncertainty"
  } else if (all(V > 0)) {
    "Tumor volume & uncertainty"
  } else {
    stop("Data must be either 0–1 (confluency) or positive (tumor volume).")
  }

  y_label_Nor <- if (all(V >= 0 & V <= 1)) {
    "Normalized confluency & uncertainty"
  } else {
    "Normalized tumor volume & uncertainty"
  }

  ## ---------- PPC plot ----------
  ppc_plot <- ggplot2::ggplot() +
    ggplot2::geom_violin(
      data = approx_data,
      ggplot2::aes(x = treatment, y = predicted_population),
      fill = "deepskyblue",
      alpha = 0.7,
      scale = "area",
      trim = FALSE
    ) +
    ggplot2::geom_point(
      data = data_v,
      ggplot2::aes(x = treat, y = V),
      colour = "red",
      size = 2
    ) +
    ggplot2::labs(
      x = "Treatment",
      y = y_label
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 2
      ),
      axis.text.x = ggplot2::element_text(
        face = "bold",
        size = 14,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(face = "bold", size = 14),
      axis.title.x = ggplot2::element_text(face = "bold", size = 16),
      axis.title.y = ggplot2::element_text(face = "bold", size = 16),
      plot.margin = ggplot2::margin(20, 20, 20, 20)
    ) +
    ggplot2::scale_y_log10()

  ## ---------- normalized PPC ----------
  ppc_plot_Nor_V <- NULL

  if (!is.null(VT0)) {

    ppc_plot_Nor_V <- ggplot2::ggplot() +
      ggplot2::geom_violin(
        data = approx_data_Nor_V,
        ggplot2::aes(x = treatment, y = predicted_population),
        fill = "deepskyblue",
        alpha = 0.7,
        scale = "area",
        trim = FALSE
      ) +
      ggplot2::geom_point(
        data = data_Nor_V,
        ggplot2::aes(x = treat, y = V),
        colour = "red",
        size = 2
      ) +
      ggplot2::labs(
        x = "Treatment",
        y = y_label_Nor
      ) +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(
          colour = "black",
          fill = NA,
          linewidth = 2
        ),
        axis.text.x = ggplot2::element_text(
          face = "bold",
          size = 14,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = ggplot2::element_text(face = "bold", size = 14),
        axis.title.x = ggplot2::element_text(face = "bold", size = 16),
        axis.title.y = ggplot2::element_text(face = "bold", size = 16),
        plot.margin = ggplot2::margin(20, 20, 20, 20)
      ) +
      ggplot2::scale_y_log10()
  }

  ## ---------- return ----------
 out <-list(
    ppc_plot_normalized = ppc_plot_Nor_V,
    ppc_plot = ppc_plot,
    ppc_data = approx_data,
    sampled_population = sampled_elements,
    type = "ppc_population"
  )
  class(out) <- "barmixR_diagnostic"
  out
}
