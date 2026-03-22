#' Rank treatments within each cell line using QTR boxplot summaries
#'
#' Visualizes quantitative treatment resistance (QTR) by ranking
#' treatments within each cell line using the summary statistics returned
#' by \code{QTRresistance()}.
#'
#' Treatments are ordered within each cell line by the posterior median
#' QTR estimate. The plot displays boxplot-style summaries using the
#' interquartile range, whiskers, and median.
#'
#' @param @param model An object of class \code{\link{barmixR_fit}} returned by \code{\link{barmixRQTR}}
#'   containing \code{names_cell}.
#'
#' @param summary_table Data frame returned by \code{QTRresistance()}.
#'   It must contain the columns:
#'   \itemize{
#'     \item \code{cell}
#'     \item \code{treat}
#'     \item \code{median}
#'     \item \code{q25}
#'     \item \code{q75}
#'     \item \code{lower_whisker}
#'     \item \code{upper_whisker}
#'   }
#'
#' @param cell_order Optional character vector giving the desired order
#'   of cell-line panels.
#'
#' @param ncol Number of columns in the combined panel layout.
#'
#' @param legend_text Optional annotation shown beneath the combined
#'   figure.
#'
#' @return An object of class \code{barmixR_result} containing:
#'   \itemize{
#'     \item \code{rank_plot}: Combined multi-panel rank plot.
#'     \item \code{rank_summary}: Summary table used for plotting,
#'       including within-cell ranks.
#'     \item \code{type}: Result type ("QTR_decision").
#'   }
#'
#' @importFrom ggplot2 ggplot aes annotate geom_vline geom_segment geom_rect geom_point
#' @importFrom ggplot2 labs theme_minimal theme
#' @importFrom ggplot2 element_rect element_line element_text margin
#' @importFrom ggplot2 coord_cartesian theme_void
#' @importFrom dplyr group_by arrange mutate ungroup row_number
#' @importFrom forcats fct_rev
#' @importFrom patchwork wrap_plots plot_layout
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
#'   sample(1:50, 20, replace = TRUE),
#'   nrow = 4
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
#' V <- c(100, 120, 300, 280)
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
#' ## ---------------------------
#' ## Fit model (fast settings)
#' ## ---------------------------
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
#'   sampled_population = ppc_pop$sampled_population
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
#' ## Decision plot
#' ## ---------------------------
#'
#' decision <- QTRDecision(
#'   model = fit,
#'   summary_table = res$summary_table
#' )
#'
#' ## visualize ranking
#' decision$rank_plot
#' 
#' @export
QTRDecision <- function(
    model,
    summary_table,
    cell_order = NULL,
    ncol = 4,
    legend_text = NULL
) {
  
  ## ---------- helpers ----------
  .stop_if <- function(cond, msg) {
    if (isTRUE(cond)) stop(msg, call. = FALSE)
  }
  
  ## ---------- validate inputs ----------
  .stop_if(is.null(model) || !inherits(model, "barmixR_fit"),
           "`model` must be an object of class 'barmixR_fit' returned by `barmixRQTR()`.")
  
  .stop_if(is.null(model$names_cell),
           "`model$names_cell` is missing.")
  
  summary_table <- as.data.frame(summary_table)
  
  required_cols <- c(
    "cell", "treat", "median", "q25", "q75",
    "lower_whisker", "upper_whisker"
  )
  missing_cols <- setdiff(required_cols, names(summary_table))
  
  .stop_if(length(missing_cols) > 0,
           paste0("Missing required columns in `summary_table`: ",
                  paste(missing_cols, collapse = ", ")))
  
  ## ---------- prepare rank summary ----------
  rank_summary <- summary_table
  rank_summary$cell <- as.character(rank_summary$cell)
  rank_summary$treat <- as.character(rank_summary$treat)
  
  rank_summary <- dplyr::group_by(rank_summary, cell)
  rank_summary <- dplyr::arrange(rank_summary, median, .by_group = TRUE)
  rank_summary <- dplyr::mutate(rank_summary, rank = dplyr::row_number())
  rank_summary <- dplyr::ungroup(rank_summary)
  
  ## ---------- per-cell plots ----------
  split_data <- split(rank_summary, rank_summary$cell)
  
  plot_list <- lapply(split_data, function(df) {
    
    df <- df[order(df$median), , drop = FALSE]
    df$treat <- factor(df$treat, levels = unique(df$treat))
    df$y <- rev(seq_len(nrow(df)))
    
    ggplot2::ggplot(df) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "black",
        linewidth = 1
      ) +
      ## whiskers
      ggplot2::geom_segment(
        ggplot2::aes(
          x = lower_whisker,
          xend = q25,
          y = y,
          yend = y
        ),
        linewidth = 0.8,
        color = "black"
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = q75,
          xend = upper_whisker,
          y = y,
          yend = y
        ),
        linewidth = 0.8,
        color = "black"
      ) +
      ## box
      ggplot2::geom_rect(
        ggplot2::aes(
          xmin = q25,
          xmax = q75,
          ymin = y - 0.25,
          ymax = y + 0.25
        ),
        fill = grDevices::adjustcolor("gray70", alpha.f = 0.2),
        color = "black",
        linewidth = 0.5
      ) +
      ## median
      ggplot2::geom_segment(
        ggplot2::aes(
          x = median,
          xend = median,
          y = y - 0.25,
          yend = y + 0.25
        ),
        linewidth = 0.9,
        color = "black"
      ) +
      ggplot2::scale_y_continuous(
        breaks = df$y,
        labels = levels(df$treat),
        expand = ggplot2::expansion(mult = c(0.05, 0.05))
      ) +
      ggplot2::labs(
        title = unique(df$cell),
        x = "Treatment resistance",
        y = "Treatment"
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.grid.major = ggplot2::element_line(color = "gray90"),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.title = ggplot2::element_text(size = 18, hjust = 0.5, color = "black", face = "bold"),
        axis.text.y = ggplot2::element_text(size = 15, color = "black"),
        axis.text.x = ggplot2::element_text(size = 15, color = "black"),
        axis.title.x = ggplot2::element_text(size = 18, color = "black"),
        axis.title.y = ggplot2::element_text(size = 18, color = "black"),
        axis.line = ggplot2::element_line(color = "black"),
        axis.ticks = ggplot2::element_line(color = "black"),
        legend.position = "none"
      )
  })
  
  ## ---------- reorder panels if requested ----------
  if (!is.null(cell_order)) {
    cell_order <- as.character(cell_order)
    keep_order <- cell_order[cell_order %in% names(plot_list)]
    plot_list <- plot_list[keep_order]
  }
  
  ## ---------- combine plots ----------
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol)
  
  if (!is.null(legend_text)) {
    legend_plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0.5, y = 0.5,
        label = legend_text,
        size = 5.5,
        hjust = 0.5,
        fontface = "plain",
        color = "black"
      ) +
      ggplot2::coord_cartesian(
        xlim = c(0, 1),
        ylim = c(0, 1),
        expand = FALSE
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.background = ggplot2::element_rect(
          fill = "gray95",
          color = "black",
          linewidth = 0.5
        ),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      )
    
    combined_plot <- (combined_plot / legend_plot) +
      patchwork::plot_layout(heights = c(1, 0.06))
  }
  
  ## ---------- return ----------
  out <- list(
    rank_plot = combined_plot,
    rank_summary = rank_summary,
    type = "QTR_decision"
  )
  class(out) <- "barmixR_result"
  out
}