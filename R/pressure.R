# Internal reference lookup; built once per R session from the bundled CSV.
.pressure_cache <- new.env(parent = emptyenv())
.reference_pressure <- function(t) {
    if (is.null(.pressure_cache$lookup)) {
      path <- system.file("extdata", "reference_pressure.csv", package = "CorPresDefl")
      d <- utils::read.csv(path)
      stopifnot(identical(names(d), c("t_ms", "P_mmHg")),
                all(is.finite(as.matrix(d))), all(diff(d$t_ms) > 0),
                all(d$P_mmHg >= 0), d$t_ms[1L] == 0,
                utils::tail(d$t_ms, 1L) == 30)
      peak <- which.max(d$P_mmHg)
      rise <- seq_len(peak)
      fall <- peak:nrow(d)
      stopifnot(peak > 1L, peak < nrow(d),
                all(diff(d$P_mmHg[rise]) >= 0),
                all(diff(d$P_mmHg[fall]) <= 0))
      up <- stats::splinefun(d$t_ms[rise], d$P_mmHg[rise], method = "monoH.FC")
      down <- stats::splinefun(d$t_ms[fall], d$P_mmHg[fall], method = "monoH.FC")
      peak_t <- d$t_ms[peak]
      .pressure_cache$lookup <- function(x) {
        y <- rep(NA_real_, length(x))
        valid <- is.finite(x) & x >= 0 & x <= 30
        before <- valid & x <= peak_t
        after <- valid & x > peak_t
        y[before] <- pmax(0, up(x[before]))
        y[after] <- pmax(0, down(x[after]))
        y
      }
    }
    .pressure_cache$lookup(t)
}
