#' Visualize quantitative treatment resistance (QTR) using bubble heatmaps
#'
#' Visualizes quantitative treatment resistance (QTR) differences between
#' cell lines and treatments using bubble heatmaps. The visualization is
#' derived from posterior summaries of treatment resistance computed with
#' \code{QTRresistance()}.
#'
#' Bubble size represents the magnitude of resistance differences, while
#' color indicates the direction of resistance change.
#'
#' @param model An object of class \code{\link{barmixR_fit}} returned by \code{\link{barmixRQTR}}
#'   containing \code{condition_count}.
#'
#' @param summary_table Data frame returned by \code{QTRresistance()}
#'   containing summary statistics for treatment resistance.
#'   The table must include the columns:
#'   \itemize{
#'     \item \code{cell}: Cell line identifier.
#'     \item \code{treat}: Treatment identifier.
#'     \item \code{mean}: Mean treatment resistance estimate.
#'     \item \code{CDF}: Cumulative distribution function value for
#'       resistance estimates.
#'   }
#'
#' @param q_up Upper quantile threshold (unused; retained for API
#'   compatibility).
#'
#' @param q_lo Lower quantile threshold (unused; retained for API
#'   compatibility).
#'
#' @param min_w Minimum bubble size used in the within-treatment heatmap.
#'
#' @param max_w Maximum bubble size used in the within-treatment heatmap.
#'
#' @param min_b Minimum bubble size used in the treatment-level heatmap.
#'
#' @param max_b Maximum bubble size used in the treatment-level heatmap.
#'
#' @return An object of class \code{barmixR_result} containing:
#'   \itemize{
#'     \item \code{heatmap_within}: Bubble heatmap showing pairwise
#'       resistance differences between cell lines within each treatment.
#'
#'     \item \code{heatmap_treatment}: Bubble heatmap summarizing
#'       resistance differences across treatments.
#'       
#'     \item \code{type}: Result type ("QTR_heatmap").
#'   }
#'
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap labs theme
#' @importFrom ggplot2 element_text element_rect margin scale_radius
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom dplyr arrange
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
#' ## Full downstream pipeline
#' ## ---------------------------
#'
#' # barcode PPC
#' ppc_bar <- ppcBarcodes(fit)
#'
#' # fraction ratios
#' frac_ratio <- fractionRatio(
#'   model = fit,
#'   sampled_fraction = ppc_bar$sampled_fraction
#' )
#'
#' # population PPC
#' ppc_pop <- ppcPopulation(fit)
#'
#' # population ratios
#' pop_ratio <- populationRatio(
#'   model = fit,
#'   sampled_elements = ppc_pop$sampled_population
#' )
#'
#' # resistance estimation
#' res <- QTRresistance(
#'   model = fit,
#'   li_sam_ratio_relative = frac_ratio$li_sam_ratio_relative,
#'   li_sam_ratio_V = pop_ratio$li_sam_ratio_V
#' )
#'
#' ## ---------------------------
#' ## Heatmap visualization
#' ## ---------------------------
#'
#' heat <- QTRheatmap(
#'   model = fit,
#'   summary_table = res$summary_table
#' )
#'
#' ## within-treatment heatmap
#' heat$heatmap_within
#'
#' ## treatment-level heatmap
#' heat$heatmap_treatment
#' 
#' @export
QTRheatmap <- function(
    model,
    summary_table,
    q_up = 0.975,
    q_lo = 0.025,
    min_w = 1,
    max_w = 8,
    min_b = 1,
    max_b = 8
) {

  ## ---------- helpers ----------
  .stop_if <- function(cond, msg) {
    if (isTRUE(cond)) stop(msg, call. = FALSE)
  }

  ## ---------- validate inputs ----------
  .stop_if(is.null(model) || !inherits(model, "barmixR_fit"),
           "`model` must be an object of class 'barmixR_fit' returned by `barmixRQTR()`.")

  .stop_if(is.null(model$condition_count),
           "`model$condition_count` is missing.")

  summary_table <- as.data.frame(summary_table)

  required_cols <- c("cell", "treat", "mean", "CDF")

 .stop_if(!all(required_cols %in% names(summary_table)), 
          "`summary_table` must contain columns: cell, treat, mean, CDF.")

  condition_count <- unique(model$condition_count)

  ## ---------- combine inputs ----------
  cdf_mean <- summary_table

  result_list <- list()

  ## ---------- compute pairwise differences ----------
  for (treatment in unique(cdf_mean$treat)) {

    subset_data <- cdf_mean[cdf_mean$treat == treatment, ]

    for (i in seq_len(nrow(subset_data))) {
      for (j in seq_len(nrow(subset_data))) {

        if (i == j) {

          direction_val <- 2 * ((1 - subset_data$CDF[i]) - 0.5)
          mean_val <- abs(subset_data$mean[i])

        } else {

          direction_val <- 0.5 * (
            (1 - subset_data$CDF[i]) -
              (1 - subset_data$CDF[j]) -
              (subset_data$CDF[i] - subset_data$CDF[j])
          )

          mean_val <- abs(subset_data$mean[i] - subset_data$mean[j])
        }

        result_list[[length(result_list) + 1]] <- data.frame(
          new = subset_data$cell[i],
          new1 = subset_data$cell[j],
          Direction = direction_val,
          treat = treatment,
          Mean = mean_val
        )
      }
    }
  }

  result1 <- do.call(rbind, result_list)

  ## ---------- factor ordering ----------
  result1$new <- factor(result1$new, levels = unique(result1$new))
  result1$new1 <- factor(result1$new1, levels = unique(result1$new1))
  result1$treat <- factor(result1$treat, levels = condition_count)

  result1 <- dplyr::arrange(result1, treat, new, new1)

  blue_side <- c("#08306B","#08519C","#2171B5","#4292C6",
               "#6BAED6","#9ECAE1","#C6DBEF","#DEEBF7")

  red_side  <- c("#FEE0D2","#FCBBA1","#FC9272","#FB6A4A",
               "#EF3B2C","#CB181D","#A50F15","#67000D")
  cols <- c(blue_side, red_side)
  ## ---------- within-treatment heatmap ----------
  bubel_heatmap_w <- ggplot2::ggplot(
    result1,
    ggplot2::aes(
      x = new,
      y = new1,
      size = Mean,
      color = Direction
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_point(shape = 21, stroke = 1.2, color = "black") +
    ggplot2::facet_wrap(~ treat) +
    ggplot2::scale_radius(range = c(min_w, max_w)) +
    ggplot2::scale_color_gradientn(
      colours = cols,
      limits = c(-1,1)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        override.aes = list(
        shape = 21,
        fill = "grey80",
        color = "black"
       )
     )
    )+
    ggplot2::labs(
      x = "Cell line",
      y = "Cell line"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold", size = 15, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(face = "bold", size = 15),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18),
      strip.text = ggplot2::element_text(face = "bold", size = 18),
      plot.margin = ggplot2::margin(20,20,40,40)
    )

  ## ---------- treatment-level heatmap ----------
  result_home <- result1[result1$new == result1$new1, ]

  bubel_heatmap_b <- ggplot2::ggplot(
    result_home,
    ggplot2::aes(
      x = treat,
      y = new1,
      size = Mean,
      color = Direction
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_point(shape = 21, stroke = 1.2, color = "black") +
    ggplot2::scale_radius(range = c(min_b, max_b)) +
    ggplot2::scale_color_gradientn(
      colours = cols,
      limits = c(-1,1)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        override.aes = list(
        shape = 21,
        fill = "grey80",
        color = "black"
       )
     )
    )+
    ggplot2::labs(
      x = "Treatment",
      y = "Cell line"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "bold", size = 15, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(face = "bold", size = 15),
      axis.title.x = ggplot2::element_text(face = "bold", size = 18),
      axis.title.y = ggplot2::element_text(face = "bold", size = 18),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 2),
      plot.margin = ggplot2::margin(20,20,40,40)
    )

  ## ---------- return ----------
  out <-list(
    heatmap_within = bubel_heatmap_w,
    heatmap_treatment = bubel_heatmap_b,
    type = "QTR_heatmap"
  )
  class(out) <- "barmixR_result"
  out
}
