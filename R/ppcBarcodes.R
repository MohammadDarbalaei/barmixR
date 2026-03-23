#' Posterior predictive checks for barcode composition
#'
#' Performs posterior predictive checks (PPC) for barcode composition
#' modeled using the Dirichlet–multinomial component of the barmixR
#' model. Observed barcode counts are compared with predicted counts
#' from the posterior distribution to assess model fit and uncertainty.
#'
#' Posterior predictive distributions are visualized using violin plots,
#' with observed counts overlaid for each cell line and treatment
#' condition.
#'
#' @param model An object of class \code{barmixR_fit} returned by \code{barmixRQTR()}
#'   containing:
#'   \itemize{
#'     \item \code{fit_dir_mult}: Fitted Dirichlet–multinomial model.
#'     \item \code{group}: Numeric identifiers for treatment groups.
#'     \item \code{condition_count}: Treatment labels for barcode data.
#'     \item \code{chains}: Number of MCMC chains used during fitting.
#'     \item \code{L}: Number of cell lines.
#'     \item \code{names_cell}: Cell line identifiers.
#'     \item \code{data_count}: Observed barcode count matrix.
#'     \item \code{K}: Number of treatment groups.
#'     \item \code{n_sam}: Number of posterior samples retained.
#'   }
#'
#' @param q_up Upper quantile threshold used for uncertainty
#'   visualization.
#'
#' @param q_lo Lower quantile threshold used for uncertainty
#'   visualization.
#'
#' @return An object of class \code{barmixR_diagnostic}, containing:
#'   \itemize{
#'     \item \code{ppc_plot}: A \code{ggplot2} posterior predictive
#'       check visualization for barcode composition.
#'
#'     \item \code{sampled_fraction}: Posterior samples of barcode
#'       compositional fractions used for downstream analysis.
#'       
#'     \item \code{type}: Result type ("ppc_barcodes"). 
#'   }
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_point facet_wrap labs theme
#' @importFrom ggplot2 element_text margin guides guide_legend scale_y_log10
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
#' ## Posterior Predictive Checks for Barcode Counts
#' ## ---------------------------
#'
#' ppc <- ppcBarcodes(fit)
#'
#' # plot PPC
#' ppc$ppc_plot
#'
#' # posterior sampled fractions
#' names(ppc$sampled_fraction)
#' 
#' 
#' @export
ppcBarcodes <- function(
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

  required_fields <- c(
    "fit_dir_mult",
    "group",
    "condition_count",
    "chains",
    "L",
    "names_cell",
    "data_count",
    "K",
    "n_sam"
  )

  missing_fields <- setdiff(required_fields, names(model))

  .stop_if(is.null(model) || !inherits(model, "barmixR_fit"),
           "`model` must be an object of class 'barmixR_fit' returned by `barmixRQTR()`.")

  .stop_if(!requireNamespace("rstan", quietly = TRUE),
           "Package `rstan` is required but not installed.")

  ## ---------- extract model components ----------
  group <- model$group
  condition_count <- as.character(model$condition_count)
  n_sam <- as.integer(model$n_sam)
  fit_dir_mult <- model$fit_dir_mult
  L <- as.integer(model$L)
  names_cell <- as.character(model$names_cell)
  y <- as.matrix(model$data_count)
  K <- as.integer(model$K)

  n <- nrow(y)
  L_cell <- length(names_cell)

  conditions <- unique(model$condition_count[order(model$group)])

  ## ---------- extract posterior samples ----------
  yrep_ab <- rstan::extract(fit_dir_mult, "y_rep")$y_rep
  softmax_fraction <- rstan::extract(fit_dir_mult, pars = "theta")$theta
  kappa <- rstan::extract(fit_dir_mult, pars = "S")$S

  ## ---------- normalize predicted fractions ----------
  for (j in seq_len(nrow(y))) {
    softmax_fraction[, j, ] <-
      t(apply(softmax_fraction[, j, ] * kappa, 1,
              function(row) row / sum(row)))
  }

  ## ---------- prepare PPC samples ----------
  final_matrix <- list()
  sampled_elements <- list()
  sampled_fraction <- list()
  fraction_matrix <- list()

  seq_list <- split(seq_len(nrow(y)), group)

  for (k in seq_len(K)) {

    unc_volume_treat <- yrep_ab[, seq_list[[k]], ]

    vectors_list <- lapply(unc_volume_treat, as.vector)

    final_matrix[[k]] <- matrix(
      unlist(vectors_list),
      nrow = dim(unc_volume_treat)[1] * length(seq_list[[k]]),
      ncol = L_cell,
      byrow = FALSE
    )

    sampled_elements[[k]] <-
      final_matrix[[k]][
        sample(seq_len(nrow(final_matrix[[k]])), n_sam),
      ]

    fraction <- softmax_fraction[, seq_list[[k]], ]
    fraction <- lapply(fraction, as.vector)

    fraction_matrix[[k]] <- matrix(
      unlist(fraction),
      nrow = dim(unc_volume_treat)[1] * length(seq_list[[k]]),
      ncol = L_cell,
      byrow = FALSE
    )
    colnames(fraction_matrix[[k]]) <- names_cell
    sampled_fraction[[k]] <-
      fraction_matrix[[k]][
        sample(seq_len(nrow(fraction_matrix[[k]])), n_sam),
      ]
  }
  names(sampled_fraction) <- conditions
  ## ---------- prepare visualization data ----------
  cell_ID <- ave(
    names_cell,
    names_cell,
    FUN = function(x) {
      if (length(x) > 1) paste0(x, seq_along(x)) else x
    }
  )

  cell_bar <- rep(
    rep(factor(cell_ID, levels = cell_ID), each = n_sam),
    length(conditions)
  )

  treatments <- rep(
    conditions,
    each = n_sam * length(cell_ID)
  )

  data_real <- data.frame(
    count = as.vector(t(y)) + 1,
    cell = rep(factor(cell_ID, levels = cell_ID), n),
    treat = rep(condition_count, each = L_cell),
    stringsAsFactors = FALSE
  )

  data_ppc <- data.frame(
    pred = unlist(sampled_elements) + 1,
    cell = cell_bar,
    treat = treatments,
    stringsAsFactors = FALSE
  )

  data_ppc$treat <- factor(data_ppc$treat, levels = conditions)

  ## ---------- plot ----------
  ppc_plot <- ggplot2::ggplot(
    data_ppc,
    ggplot2::aes(x = cell, y = pred, fill = cell)
  ) +
    ggplot2::geom_violin(
      scale = "width",
      trim = FALSE
    ) +
    ggplot2::facet_wrap(~ treat) +
    ggplot2::geom_point(
      data = data_real,
      ggplot2::aes(x = cell, y = count)
    ) +
    ggplot2::labs(
      x = "Cell",
      y = "Absolute abundance & uncertainty",
      fill = "Cell"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        face = "bold",
        size = 18,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = ggplot2::element_text(face = "bold", size = 18),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18),
      strip.text = ggplot2::element_text(face = "bold", size = 18),
      plot.margin = ggplot2::margin(20, 20, 40, 40),
      legend.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 15)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::scale_y_log10()

  ## ---------- return ----------
  out <-list(
    ppc_plot = ppc_plot,
    sampled_fraction = sampled_fraction,
    type = "ppc_barcodes"
  )
  class(out) <- "barmixR_diagnostic"
  out
}
