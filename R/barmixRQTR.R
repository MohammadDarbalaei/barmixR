#' Fit the barmixR model for quantitative treatment resistance (QTR)
#'
#' Fits the Bayesian barmixR framework for joint modeling of barcode
#' composition and population size measurements in pooled cell
#' populations. The model integrates barcode sequencing counts with
#' population size data to estimate quantitative treatment resistance
#' (QTR) for individual cell lines.
#'
#' Population size measurements correspond to tumor volume in
#' \emph{in vivo} experiments or cellular confluency in \emph{in vitro}
#' assays. Barcode counts are modeled using a Dirichlet–multinomial
#' distribution to account for compositional sequencing data, while
#' population size measurements are modeled using a log-normal
#' likelihood for tumor volume (\emph{in vivo}) or a beta likelihood
#' for confluency (\emph{in vitro}).
#'
#' Posterior inference is performed using Hamiltonian Monte Carlo
#' implemented in \pkg{rstan}.
#'
#' @param data A list containing the following components:
#'   \itemize{
#'     \item \code{data_count}: A matrix or data frame of barcode read
#'       counts. Rows correspond to experimental replicates
#'       (e.g., \code{DMSO_1}) and columns correspond to genotypes or
#'       mutations (e.g., \code{PTEN}, \code{TSC2}).
#'
#'     \item \code{condition_count}: A factor or character vector
#'       specifying treatment groups corresponding to rows of
#'       \code{data_count}.
#'
#'     \item \code{V}: A numeric vector of population size measurements,
#'       representing tumor volume in \emph{in vivo} experiments or
#'       confluency in \emph{in vitro} assays.
#'
#'     \item \code{condition_v}: A factor or character vector specifying
#'       treatment groups corresponding to entries of \code{V}.
#'
#'     \item \code{cell_line}: A factor or character vector of cell line
#'       identifiers corresponding to the columns of \code{data_count}.
#'
#'     \item \code{VT0}: Optional numeric vector of baseline population
#'       size values used for normalization (tumor volume in
#'       \emph{in vivo} experiments or confluency in \emph{in vitro}
#'       assays).
#'
#'     \item \code{pre_composition}: Optional matrix or data frame used
#'       to adjust barcode compositions prior to model fitting.
#'   }
#'
#' @param control_group Optional label specifying the control treatment
#'   group used for comparative analysis.
#'
#' @param dispersion A list controlling dispersion priors in the
#'   Dirichlet–multinomial model:
#'   \itemize{
#'     \item \code{psi_mean}: Mean for the log-normal prior on the
#'       precision parameter \code{psi} (default: 0).
#'
#'     \item \code{psi_sd}: Standard deviation for the log-normal prior
#'       on \code{psi} (default: 0.5).
#'
#'     \item \code{varphi_mean}: Mean for the log-normal prior on the
#'       variance parameter \code{varphi} (default: 0).
#'
#'     \item \code{varphi_sd}: Standard deviation for the log-normal
#'       prior on \code{varphi} (default: 0.01).
#'   }
#'
#' @param control A list of control parameters:
#'   \itemize{
#'     \item \code{chains}: Number of MCMC chains (default: 3).
#'
#'     \item \code{iter_count}: Iterations used to fit the barcode
#'       composition model (default: 1000).
#'
#'     \item \code{iter_V}: Iterations used to fit the population size
#'       model (default: 10000).
#'
#'     \item \code{cores}: Number of CPU cores used for Stan sampling
#'       (default: 3).
#'
#'     \item \code{warmup}: Optional number of warm-up iterations passed
#'       to Stan. If \code{NULL}, Stan default behavior is used.
#'   }
#'
#' @param time_d Duration of the experiment used for estimating
#'   treatment resistance.
#'
#' @param q_up Upper quantile threshold used in downstream filtering
#'   (default: 0.975).
#'
#' @param q_lo Lower quantile threshold used in downstream filtering
#'   (default: 0.025).
#'
#' @param n_sam Number of posterior samples retained for downstream
#'   analysis (defaults to half of
#'   \code{min(iter_count, iter_V)}).
#'
#' @param in_vivo Logical indicating the experimental setting.
#'   If \code{TRUE}, population size corresponds to tumor volume in
#'   \emph{in vivo} experiments and is modeled using a log-normal
#'   likelihood. If \code{FALSE}, population size corresponds to
#'   confluency in \emph{in vitro} assays and is modeled using a beta
#'   likelihood.
#'
#' @param seed Optional integer seed used for reproducible Stan sampling.
#'
#' @param verbose Logical; if \code{TRUE}, informational messages are
#'   printed during model fitting.
#'
#' @param BPPARAM Optional \code{\link[BiocParallel]{BiocParallelParam}}
#'   object used to determine parallel workers for Stan sampling.
#'
#' @param ... Reserved for future extensions.
#'
#' @return An object of class \code{barmixR_fit}. This is a list containing:
#'   \itemize{
#'     \item \code{fit_dir_mult}: Fitted Dirichlet–multinomial model
#'       for barcode counts.
#'
#'     \item \code{fit_V}: Fitted model for population size
#'       (tumor volume in \emph{in vivo} experiments or confluency in
#'       \emph{in vitro} assays).
#'
#'     \item \code{L}: Number of unique cell lines.
#'
#'     \item \code{iter_count}, \code{iter_V}, \code{chains}: MCMC
#'       parameters used for model fitting.
#'
#'     \item \code{data_count}, \code{V}, \code{VT0}: Input data used
#'       after preprocessing and filtering.
#'
#'     \item \code{group}, \code{group_V}: Numeric identifiers for
#'       treatment groups.
#'
#'     \item \code{condition_count}, \code{condition_v}: Treatment
#'       labels for barcode counts and population size data.
#'
#'     \item \code{control_group}: Identifier for the control group.
#'
#'     \item \code{n_sam}: Number of posterior samples retained for
#'       downstream analysis.
#'
#'     \item \code{in_vivo}: Logical indicating whether the model was
#'       fitted using the \emph{in vivo} likelihood.
#'   }
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
#' ## inspect result
#' names(fit)
#' 
#' @export
barmixRQTR <- function(
    data,
    control_group = NULL,
    dispersion = list(psi_mean = 0, psi_sd = 0.5, varphi_mean = 0, varphi_sd = 0.01),
    control = list(chains = 3, iter_count = 1000, iter_V = 10000, cores = 3, warmup = NULL),
    time_d,
    q_up = 0.975,
    q_lo = 0.025,
    n_sam = NULL,
    in_vivo = TRUE,
    seed = NULL,
    verbose = FALSE,
    BPPARAM = NULL,
    ...
) {
  ## ---------- helpers ----------
  .stop_if <- function(cond, msg) {
    if (isTRUE(cond)) stop(msg, call. = FALSE)
  }
  
  .as_chr <- function(x) {
    if (is.factor(x)) as.character(x) else x
  }
  
  ## ---------- required args ----------
  .stop_if(missing(time_d) || length(time_d) != 1 || !is.finite(time_d),
           "`time_d` is required and must be a finite numeric scalar.")
  
  .stop_if(is.null(data) || !is.list(data), "`data` must be a list.")
  
  required_fields <- c("data_count", "condition_count", "V", "condition_v", "cell_line")
  missing_fields <- setdiff(required_fields, names(data))
  .stop_if(length(missing_fields) > 0,
           paste0("Missing required fields in `data`: ", paste(missing_fields, collapse = ", ")))
  
  ## ---------- coerce & validate core inputs ----------
  data_count <- data$data_count
  .stop_if(!(is.matrix(data_count) || is.data.frame(data_count)),
           "`data$data_count` must be a matrix or data.frame.")
  data_count <- as.matrix(data_count)
  .stop_if(anyNA(data_count), "`data$data_count` contains NA values.")
  .stop_if(any(data_count < 0), "`data$data_count` must be non-negative.")
  storage.mode(data_count) <- "numeric"
  
  condition_count <- .as_chr(data$condition_count)
  condition_v <- .as_chr(data$condition_v)
  cell_line <- .as_chr(data$cell_line)
  
  .stop_if(length(condition_count) != nrow(data_count),
           "`data$condition_count` must have the same length as nrow(data$data_count).")
  
  .stop_if(length(cell_line) != ncol(data_count),
           "`data$cell_line` must have the same length as ncol(data$data_count).")
  
  V <- as.numeric(data$V)
  .stop_if(length(condition_v) != length(V),
           "`data$condition_v` must have the same length as `data$V`.")
  .stop_if(anyNA(V) || anyNA(condition_v),
           "`data$V` and `data$condition_v` must not contain NA values.")
  
  VT0 <- NULL
  if (!is.null(data$VT0)) {
    VT0 <- as.numeric(data$VT0)
    .stop_if(anyNA(VT0), "`data$VT0` contains NA values.")
  }
  
  ## ---------- control validation ----------
  control <- modifyList(
    list(chains = 3L, iter_count = 1000L, iter_V = 10000L, cores = 3L, warmup = NULL),
    control
  )
  
  .stop_if(length(control$chains) != 1 || !is.finite(control$chains) || as.integer(control$chains) <= 0,
           "`control$chains` must be a positive integer.")
  control$chains <- as.integer(control$chains)
  
  .stop_if(length(control$iter_count) != 1 || !is.finite(control$iter_count) || as.integer(control$iter_count) <= 0,
           "`control$iter_count` must be a positive integer.")
  control$iter_count <- as.integer(control$iter_count)
  
  .stop_if(length(control$iter_V) != 1 || !is.finite(control$iter_V) || as.integer(control$iter_V) <= 0,
           "`control$iter_V` must be a positive integer.")
  control$iter_V <- as.integer(control$iter_V)
  
  if (!is.null(BPPARAM)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("`BPPARAM` was provided but BiocParallel is not installed.", call. = FALSE)
    }
    control$cores <- as.integer(max(1L, BiocParallel::bpnworkers(BPPARAM)))
  } else {
    .stop_if(length(control$cores) != 1 || !is.finite(control$cores) || as.integer(control$cores) <= 0,
             "`control$cores` must be a positive integer.")
    control$cores <- as.integer(control$cores)
  }
  
  if (!is.null(control$warmup)) {
    .stop_if(length(control$warmup) != 1 || !is.finite(control$warmup) ||
               as.integer(control$warmup) <= 0,
             "`control$warmup` must be a positive integer.")
    .stop_if(as.integer(control$warmup) >= min(control$iter_count, control$iter_V),
             "`control$warmup` must be smaller than both `iter_count` and `iter_V`.")
    control$warmup <- as.integer(control$warmup)
  }
  
  if (!is.null(seed)) {
    .stop_if(length(seed) != 1 || !is.finite(seed), "`seed` must be a single finite number.")
    seed <- as.integer(seed)
  }
  
  ## ---------- n_sam ----------
  if (is.null(n_sam)) {
    n_sam <- floor(min(control$iter_count, control$iter_V) / 2)
  } else {
    .stop_if(length(n_sam) != 1 || !is.finite(n_sam) || as.integer(n_sam) <= 0,
             "`n_sam` must be a positive integer.")
    n_sam <- as.integer(min(as.integer(n_sam), min(control$iter_count, control$iter_V)))
  }
  
  ## ---------- dispersion ----------
  dispersion <- modifyList(
    list(psi_mean = 0, psi_sd = 0.5, varphi_mean = 0, varphi_sd = 0.01),
    dispersion
  )
  .stop_if(any(!is.finite(unlist(dispersion))), "`dispersion` values must be finite numerics.")
  
  dispersion$psi_mean <- as.numeric(dispersion$psi_mean)
  dispersion$psi_sd <- as.numeric(dispersion$psi_sd)
  dispersion$varphi_mean <- as.numeric(dispersion$varphi_mean)
  dispersion$varphi_sd <- as.numeric(dispersion$varphi_sd)
  
  .stop_if(dispersion$psi_sd <= 0, "`dispersion$psi_sd` must be > 0.")
  .stop_if(dispersion$varphi_sd <= 0, "`dispersion$varphi_sd` must be > 0.")
  
  ## ---------- optional pre_composition adjustment ----------
  if (!is.null(data$pre_composition)) {
    pre_comp <- data$pre_composition
    .stop_if(!(is.matrix(pre_comp) || is.data.frame(pre_comp)),
             "`data$pre_composition` must be a matrix or data.frame when provided.")
    pre_comp <- as.matrix(pre_comp)
    .stop_if(anyNA(pre_comp), "`data$pre_composition` contains NA values.")
    
    if (!all(dim(pre_comp) == dim(data_count))) {
      .stop_if(is.null(rownames(data_count)) || is.null(rownames(pre_comp)),
               "`data$pre_composition` must have the same dimensions as `data$data_count`, or matching row names.")
      .stop_if(!all(rownames(data_count) %in% rownames(pre_comp)),
               "All row names of `data$data_count` must be present in `data$pre_composition`.")
      pre_comp <- pre_comp[rownames(data_count), , drop = FALSE]
      .stop_if(!all(dim(pre_comp) == dim(data_count)),
               "`data$pre_composition` could not be aligned to `data$data_count`.")
    }
    
    adj <- data_count
    for (i in seq_len(nrow(adj))) {
      row_sum <- sum(adj[i, ])
      .stop_if(!is.finite(row_sum) || row_sum <= 0,
               "Each row of `data$data_count` must have a positive total count.")
      frac <- adj[i, ] / row_sum
      adj_i <- round(adj[i, ] * frac / pre_comp[i, ])
      adj_i[!is.finite(adj_i)] <- 1
      adj[i, ] <- pmax(1, adj_i)
    }
    data_count <- adj
  } else {
    if (isTRUE(verbose)) {
      message("Each barcode distributed identically among conditions.")
    }
  }
  
  ## ---------- grouping ----------
  factor_group <- factor(condition_count, levels = unique(condition_count))
  group <- as.integer(factor_group)
  K <- nlevels(factor_group)
  L <- length(unique(cell_line))
  n <- nrow(data_count)
  
  ## ---------- control group ----------
  if (!is.null(control_group)) {
    .stop_if(length(control_group) != 1, "`control_group` must be a single value.")
    control_group <- as.character(control_group)
    .stop_if(!(control_group %in% unique(condition_count)),
             "`control_group` must be present in `data$condition_count`.")
  }
  
  ## ---------- Stan availability ----------
  .stop_if(!requireNamespace("rstan", quietly = TRUE),
           "Package `rstan` is required but not installed.")
  
  ## ---------- data for count model ----------
  dlist <- list(
    y = data_count,
    group = group,
    K = K,
    C = length(cell_line),
    n = n,
    psi_mean = dispersion$psi_mean,
    psi_sd = dispersion$psi_sd,
    varphi_mean = dispersion$varphi_mean,
    varphi_sd = dispersion$varphi_sd
  )
  
  stan_seed <- seed
  if (is.null(stan_seed)) {
    stan_seed <- as.integer(sample.int(.Machine$integer.max, 1))
    if (isTRUE(verbose)) {
      message("No `seed` provided; using a randomly generated Stan seed.")
    }
  }
  
  ## ---------- fit Dirichlet-multinomial model ----------
  stan_file_counts <- system.file("stan/barcode.stan", package = "barmixR")
  
  .stop_if(!nzchar(stan_file_counts),
           "Stan file for the count model was not found in package `barmixR`.")
  
  stan_args <- list(
    file = stan_file_counts,
    data = dlist,
    chains = control$chains,
    iter = control$iter_count,
    cores = control$cores,
    seed = stan_seed,
    control = list(adapt_delta = 0.999, stepsize = 0.5, max_treedepth = 18)
  )
  
  if (!is.null(control$warmup)) {
    stan_args$warmup <- control$warmup
  }
  
  fit <- do.call(rstan::stan, stan_args)
  
  ## ---------- treatment consistency ----------
  .stop_if(!all(unique(condition_count) %in% unique(condition_v)),
           "Error: Treatments in `condition_v` must include all treatments in `condition_count`.")
  
  ## ---------- prepare V data ----------
  keep <- condition_v %in% unique(condition_count)
  V_use <- as.numeric(V[keep])
  cond_v_use <- condition_v[keep]
  
  .stop_if(length(V_use) == 0,
           "No tumor volume/confluency entries match the treatments in `condition_count`.")
  
  group_V_fac <- factor(cond_v_use, levels = levels(factor_group))
  group_V <- as.integer(group_V_fac)
  
  dlist_volume <- list(
    V = V_use,
    group = group_V,
    K = nlevels(group_V_fac),
    n = length(group_V)
  )
  
  ## ---------- fit volume/confluency model ----------
  if (isTRUE(in_vivo)) {
    .stop_if(any(V_use <= 0),
             "Error: Tumor volume values must be strictly positive for the log-normal likelihood.")
    
    stan_file_V <- system.file("stan/tumor_volume.stan", package = "barmixR")
    .stop_if(!nzchar(stan_file_V),
             "Stan file for the in vivo volume model was not found in package `barmixR`.")
    
    stan_args_V <- list(
      file = stan_file_V,
      data = dlist_volume,
      chains = control$chains,
      iter = control$iter_V,
      cores = control$cores,
      seed = stan_seed,
      control = list(max_treedepth = 18)
    )
    
    if (!is.null(control$warmup)) {
      stan_args_V$warmup <- control$warmup
    }
    
    fit_V <- do.call(rstan::stan, stan_args_V)
  } else {
    .stop_if(any(V_use <= 0 | V_use >= 1),
             "Error: Confluency values must be between 0 and 1 for the beta likelihood.")
    
    stan_file_V <- system.file("stan/confluency.stan", package = "barmixR")
    .stop_if(!nzchar(stan_file_V),
             "Stan file for the in vitro confluency model was not found in package `barmixR`.")
    
    stan_args_V <- list(
      file = stan_file_V,
      data = dlist_volume,
      chains = control$chains,
      iter = control$iter_V,
      cores = control$cores,
      seed = stan_seed,
      control = list(adapt_delta = 0.999, stepsize = 0.5, max_treedepth = 12)
    )
    
    if (!is.null(control$warmup)) {
      stan_args_V$warmup <- control$warmup
    }
    
    fit_V <- do.call(rstan::stan, stan_args_V)
  }
  
  ## ---------- return ----------
  out <-list(
    fit_dir_mult = fit,
    fit_V = fit_V,
    L = L,
    iter_count = control$iter_count,
    iter_V = control$iter_V,
    chains = control$chains,
    data_count = data_count,
    V = V_use,
    VT0 = VT0,
    K = K,
    time_d = time_d,
    names_cell = cell_line,
    group = group,
    group_V = group_V,
    condition_count = condition_count,
    condition_v = condition_v,
    control_group = control_group,
    n_sam = n_sam,
    in_vivo = in_vivo
  )
  class(out) <- "barmixR_fit"
  out
}
