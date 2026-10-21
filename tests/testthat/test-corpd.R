fixture <- function() data.frame(
  A1t = c(7.5, 8), HCt = c(17, 17.5), A2t = c(22, 22.5),
  bIOP = c(14, 15), A1DeflVol = c(.2, .25),
  HCDeflVol = c(4, 4.5), A2DeflVol = c(.4, .5))

test_that("bundled waveform retains exact knots and support", {
  d <- read.csv(system.file("extdata", "reference_pressure.csv", package = "CorPresDefl"))
  expect_equal(nrow(d), 601L)
  expect_equal(.reference_pressure(d$t_ms), d$P_mmHg, tolerance = 1e-12)
  expect_equal(.reference_pressure(c(0, 18.2, 30)), c(0, 94.579257416, 0))
  expect_true(all(is.na(.reference_pressure(c(NA_real_, Inf, -1, 31)))))
  expect_length(.reference_pressure(numeric()), 0)
  tt <- seq(0, 30, length.out = 6001)
  p <- .reference_pressure(tt)
  expect_true(all(diff(p[tt <= 18.2]) >= -1e-12))
  expect_true(all(diff(p[tt >= 18.2]) <= 1e-12))
})

test_that("three endpoint slopes agree with independent calculation", {
  p <- .reference_pressure(c(7.5, 17, 22)) - 14
  expected <- c((exp(p[1]/4)-exp(-14/4))/.2,
    log(p[2]/p[1])/3.8, (exp(p[3]/20)-exp(p[2]/20))/(.4-4))
  z <- corpd(7.5, 17, 22, 14, .2, 4, .4)
  expect_equal(unname(unlist(z[1:3])), expected, tolerance = 1e-12)
  expect_identical(z$status, "OK")
  expect_equal(attr(z, "pressure_reference")$c1_mmHg, 4)
  expect_identical(attr(z, "pressure_reference")$id,
                   "corpresdefl_reference_pressure_v1")
  expect_identical(attr(z, "pressure_reference")$doi,
                   "10.1080/10255842.2015.1077515")
  expect_identical(names(z), c("k_initial_a1", "k_a1_hc", "k_hc_a2", "status"))
  expect_identical(unname(attr(z, "parameter_units")), rep("1/(mm^3)", 3))
})

test_that("volume units are explicit and inverse-volume scaling is correct", {
  a <- corpd(7.5, 17, 22, 14, .2, 4, .4)
  expect_equal(a, corpd(7.5, 17, 22, 14, .2, 4, .4, measure = "Vol"))
  for (m in c("Amp", "Area", "V")) {
    expect_error(corpd(7.5, 17, 22, 14, .2, 4, .4, measure = m), "Vol")
    expect_error(fixture() |> corpd(measure = m), "Vol")
  }
  b <- corpd(7.5, 17, 22, 14, 2, 40, 4)
  expect_equal(as.matrix(b[1:3]), as.matrix(a[1:3])/10)
})

test_that("native pipe appends outputs, preserving rows, columns and IDs", {
  d <- fixture(); d$id <- c("s1", "s2")
  z <- d |> corpd()
  expect_equal(z[names(d)], d, ignore_attr = TRUE)
  expect_true(all(z$Vol_status == "OK"))
  v <- corpd(d$A1t, d$HCt, d$A2t, d$bIOP, d$A1DeflVol, d$HCDeflVol, d$A2DeflVol)
  expect_equal(z$Vol_k_a1_hc, v$k_a1_hc)
  both <- d |> corpd() |> corpd(prefix = "repeat_")
  expect_equal(both$repeat_k_a1_hc, both$Vol_k_a1_hc)
  expect_length(attr(both, "corpd_calculations"), 2)
  expect_error(z |> corpd(), "already exist")
  expect_silent(z |> corpd(prefix = "repeat_"))
})

test_that("custom columns and zero rows work", {
  d <- fixture(); names(d)[5:7] <- c("v1", "vh", "v2")
  expect_silent(d |> corpd(measure = "Vol", columns = c(x_a1 = "v1", x_hc = "vh", x_a2 = "v2")))
  z <- fixture()[FALSE, ] |> corpd()
  expect_equal(nrow(z), 0L)
  expect_equal(ncol(z), 11L)
  expect_error(fixture()[1:4] |> corpd(), "Missing columns")
  expect_error(fixture() |> corpd(columns = c(unknown = "A1t")), "mapping")
})

test_that("only scalars recycle; wrong types and unknown arguments error", {
  z <- corpd(c(7.5, 8), 17, 22, 14, .2, 4, .4)
  expect_equal(nrow(z), 2L)
  expect_error(corpd(c(7, 8), c(17, 18, 19), 22, 14, .1, .95, .15), "length")
  expect_error(corpd("7.5", 17, 22, 14, .1, .95, .15), "numeric")
  expect_error(corpd(7.5, 17, 22, 14, .1, .95, .15, typo = 1), "Unrecognized")
  expect_error(corpd(7.5, 17, 22, 14, .1, .95, .15, measure = "DeflAmp"))
})

test_that("invalid rows are flagged, not silently dropped", {
  d <- fixture()[rep(1, 7), ]
  d$A1t[2] <- NA; d$A2t[3] <- 31; d$A1t[4] <- 20
  d$bIOP[5] <- 95; d$HCDeflVol[6] <- .1; d$bIOP[7] <- 0
  expect_warning(z <- d |> corpd(), "6 of 7")
  expect_equal(nrow(z), 7)
  expect_identical(z$Vol_status[1], "OK")
  expect_true(all(is.na(z[2:7, paste0("Vol_", c("k_initial_a1", "k_a1_hc", "k_hc_a2"))])))
  expect_match(z$Vol_status[2], "non_finite")
  expect_match(z$Vol_status[5], "non_positive_log")
  expect_error(d |> corpd(invalid = "error"), "6 of 7")
  expect_silent(d |> corpd(invalid = "quiet"))
})
