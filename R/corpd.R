#' Three phase-specific net pressure-volume slopes
#'
#' Compute transformed endpoint slopes for Initial-A1, A1-HC,
#' and HC-A2 using a fixed reference pressure waveform and the bIOP from
#' the same examination. Vector and data-frame inputs are supported.
#'
#' @param A1t First-applanation time in ms, or a data frame for the data-frame
#'   method. Data frames are suitable as the left-hand side of `|>`.
#' @param HCt,A2t Highest-concavity and second-applanation times in ms.
#' @param bIOP Positive biomechanically corrected IOP in mmHg from the same
#'   examination.
#' @param x_a1,x_hc,x_a2 Previously calculated DeflVol at the three events,
#'   in cubic millimetres (microlitres).
#' @param measure The input measure, `"Vol"`.
#' @param invalid Row-validation handling: `"warn"` (default), `"error"`, or
#'   `"quiet"`. Invalid rows retain their position with three NA results and a
#'   semicolon-separated `status`. Invalid argument types/lengths always error.
#' @param columns For data frames, an optional named character vector mapping
#'   any of `A1t`, `HCt`, `A2t`, `bIOP`, `x_a1`, `x_hc`, `x_a2` to column names.
#'   Defaults use `A1t`, `HCt`, `A2t`, `bIOP`, `A1DeflVol`,
#'   `HCDeflVol`, and `A2DeflVol`.
#' @param prefix Prefix for result columns appended to a data frame. Defaults
#'   to `"Vol_"`. Result column names must be distinct from existing names.
#' @param ... Arguments passed to methods. Unrecognized arguments error.
#'
#' @details
#' All seven numeric inputs must have a common length or length one; only
#' scalars are recycled. Inputs must be finite with
#' `0 < A1t < HCt < A2t <= 30`, positive bIOP and positive deflections,
#' and `x_hc > x_a1` and `x_hc > x_a2`. Both A1 and HC net pressures must
#' be positive for logarithms. Zero-row data frames return zero-row results.
#'
#' Reference pressure uses split Fritsch-Carlson monotone cubic Hermite
#' interpolation. Net pressure is reference pressure minus bIOP.
#' Initial deflection is zero and initial net pressure is minus bIOP.
#' The three pressure transformations are `exp(P_net / 4)`,
#' `log(P_net / (1 mmHg))`, and `exp(P_net / 20)`.
#' Slopes are signed endpoint differences of transformed pressure divided by
#' signed volume differences. Slopes have units of inverse cubic millimetres.
#'
#' @return Numeric/vector method: a data frame with `k_initial_a1`, `k_a1_hc`,
#'   `k_hc_a2`, and `status`
#'   (`"OK"` or reasons). Data-frame method: the original rows/columns with
#'   prefixed results appended. Vector-result attributes `measure`,
#'   `parameter_units`, and `pressure_reference` describe the calculation;
#'   data frames retain metadata in `corpd_calculations`.
#' @references
#' Joda AA, Shervin MMS, Kook D, Elsheikh A (2016).
#' Development and validation of a correction equation for Corvis tonometry.
#' Computer Methods in Biomechanics and Biomedical Engineering 19(9):943-953.
#' \doi{10.1080/10255842.2015.1077515}.
#'
#' Ren Q, Deng Z, Chu Z, Pan Z, Wan P (2026).
#' DeflVol: a model-derived volumetric descriptor of corneal deflection at
#' highest concavity under Corvis ST air-puff loading.
#' Frontiers in Bioengineering and Biotechnology 14:1933733.
#' \doi{10.3389/fbioe.2026.1933733}.
#' @examples
#' corpd(7.5, 17, 22, 14, 0.2, 4, 0.4)
#' corpd(c(7.5, 8), c(17, 17.5), c(22, 22.5), c(14, 15),
#'       c(.2, .25), c(4, 4.5), c(.4, .5))
#' df <- data.frame(A1t = 7.5, HCt = 17, A2t = 22, bIOP = 14,
#'                  A1DeflVol = .2, HCDeflVol = 4, A2DeflVol = .4)
#' df |> corpd()
#' @export
corpd <- function(A1t, ...) UseMethod("corpd")

#' @rdname corpd
#' @export
corpd.default <- function(A1t, HCt, A2t, bIOP, x_a1, x_hc, x_a2,
                          measure = "Vol",
                          invalid = c("warn", "error", "quiet"), ...) {
  if (length(list(...))) stop("Unrecognized arguments in `...`.", call. = FALSE)
  .check_volume_measure(measure)
  invalid <- match.arg(invalid)
  inputs <- list(A1t = A1t, HCt = HCt, A2t = A2t, bIOP = bIOP,
                 x_a1 = x_a1, x_hc = x_hc, x_a2 = x_a2)
  numeric_vector <- function(x) is.numeric(x) && is.null(dim(x)) && !is.object(x)
  if (!all(vapply(inputs, numeric_vector, logical(1)))) {
    stop("All seven inputs must be plain numeric vectors.", call. = FALSE)
  }
  sizes <- lengths(inputs)
  n <- if (any(sizes == 0L)) 0L else max(sizes)
  if (!all(sizes %in% c(1L, n))) {
    stop("Input lengths must agree; only length-one scalars are recycled.", call. = FALSE)
  }
  v <- lapply(inputs, rep_len, length.out = n)
  status <- rep("", n)
  flag <- function(bad, reason) {
    at <- which(bad %in% TRUE)
    status[at] <<- ifelse(nzchar(status[at]),
                         paste(status[at], reason, sep = ";"), reason)
  }
  finite <- Reduce(`&`, lapply(v, is.finite))
  flag(!finite, "non_finite_input")
  flag(finite & (v$A1t <= 0 | v$HCt > 30 | v$A2t > 30 |
                 v$A1t >= v$HCt | v$HCt >= v$A2t), "invalid_event_times")
  flag(finite & v$bIOP <= 0, "non_positive_bIOP")
  flag(finite & (v$x_a1 <= 0 | v$x_hc <= 0 | v$x_a2 <= 0 |
                 v$x_hc <= v$x_a1 | v$x_hc <= v$x_a2), "invalid_deflection_order")
  p_a <- .reference_pressure(v$A1t) - v$bIOP
  p_h <- .reference_pressure(v$HCt) - v$bIOP
  p_b <- .reference_pressure(v$A2t) - v$bIOP
  flag(finite & (p_a <= 0 | p_h <= 0), "non_positive_log_pressure")
  names_out <- c("k_initial_a1", "k_a1_hc", "k_hc_a2")
  out <- as.data.frame(matrix(NA_real_, nrow = n, ncol = 3L))
  names(out) <- names_out
  ok <- which(!nzchar(status))
  if (length(ok)) {
    xa <- v$x_a1[ok]; xh <- v$x_hc[ok]; xb <- v$x_a2[ok]
    y0 <- exp(-v$bIOP[ok] / 4)
    ya <- exp(p_a[ok] / 4)
    la <- log(p_a[ok]); lh <- log(p_h[ok])
    eh <- exp(p_h[ok] / 20); eb <- exp(p_b[ok] / 20)
    out[ok, ] <- cbind((ya - y0) / xa, (lh - la) / (xh - xa),
                       (eb - eh) / (xb - xh))
    bad <- rep(FALSE, n)
    bad[ok] <- !apply(is.finite(as.matrix(out[ok, , drop = FALSE])), 1L, all)
    flag(bad, "non_finite_result")
  }
  bad <- nzchar(status)
  out[bad, ] <- NA_real_
  out$status <- status
  out$status[!bad] <- "OK"
  if (any(bad)) {
    msg <- sprintf("%d of %d rows invalid; inspect `status`. No rows were removed.", sum(bad), n)
    if (invalid == "error") stop(msg, call. = FALSE)
    if (invalid == "warn") warning(msg, call. = FALSE)
  }
  attr(out, "measure") <- measure
  attr(out, "parameter_units") <- stats::setNames(rep("1/(mm^3)", 3), names_out)
  attr(out, "pressure_reference") <- list(id = "corpresdefl_reference_pressure_v1",
    doi = "10.1080/10255842.2015.1077515",
    csv_sha256 = "9ceedbcfee7dc1ee37424b78a3b80b99700bd7b7cd6499c5ea287f2619043f2f",
    interpolation = "split monoH.FC", c1_mmHg = 4, c2_mmHg = 20,
    log_reference_mmHg = 1, initial_deflection = 0)
  out
}

#' @rdname corpd
#' @export
corpd.data.frame <- function(A1t, measure = "Vol",
                             columns = NULL, prefix = NULL,
                             invalid = c("warn", "error", "quiet"), ...) {
  if (length(list(...))) stop("Unrecognized arguments in `...`.", call. = FALSE)
  .check_volume_measure(measure)
  invalid <- match.arg(invalid)
  if (anyDuplicated(names(A1t))) stop("Data-frame column names must be unique.", call. = FALSE)
  map <- c(A1t = "A1t", HCt = "HCt", A2t = "A2t", bIOP = "bIOP",
           x_a1 = paste0("A1Defl", measure), x_hc = paste0("HCDefl", measure),
           x_a2 = paste0("A2Defl", measure))
  if (!is.null(columns)) {
    if (!is.character(columns) || is.null(names(columns)) ||
        anyNA(columns) || any(!nzchar(columns)) || anyNA(names(columns)) ||
        anyDuplicated(names(columns)) || !all(names(columns) %in% names(map))) {
      stop("`columns` must be a named character mapping of known input names.", call. = FALSE)
    }
    map[names(columns)] <- columns
  }
  missing <- setdiff(unname(map), names(A1t))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (is.null(prefix)) prefix <- paste0(measure, "_")
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) {
    stop("`prefix` must be a single character string.", call. = FALSE)
  }
  args <- lapply(map, function(nm) A1t[[nm]])
  result <- do.call(corpd.default, c(args, list(measure = measure, invalid = invalid)))
  new_names <- paste0(prefix, names(result))
  if (any(new_names %in% names(A1t))) {
    stop("Result columns already exist; choose a different `prefix`.", call. = FALSE)
  }
  out <- as.data.frame(A1t)
  out[new_names] <- result
  meta <- list(measure = attr(result, "measure"), units = attr(result, "parameter_units"),
               pressure_reference = attr(result, "pressure_reference"))
  attr(out, "corpd_calculations") <- c(attr(A1t, "corpd_calculations"), stats::setNames(list(meta), prefix))
  out
}

.check_volume_measure <- function(measure) {
  if (!identical(measure, "Vol")) {
    stop('`measure` must be "Vol"; supply DeflVol in mm^3 at A1, HC, and A2.',
         call. = FALSE)
  }
}
