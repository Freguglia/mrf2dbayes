library(mrf2d)

# Shared fixtures --------------------------------------------------------
# Small field to keep tests fast
set.seed(42)
z_small <- matrix(sample(0:2, 20 * 20, replace = TRUE), 20, 20)
R <- mrfi(1)
fam <- "onepar"

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

test_that("initialize returns an R6 object of class MRFPseudoBayes", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  expect_s3_class(s, "MRFPseudoBayes")
  expect_true(R6::is.R6(s))
})

test_that("initialize: n_samples starts at 0 and samples is NULL", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  expect_equal(s$n_samples, 0L)
  expect_null(s$samples)
})

test_that("initialize: init_theta = 'zero' sets theta to all-zeros", {
  s <- MRFPseudoBayes$new(z_small, R, fam, init_theta = "zero")
  expect_true(all(s$theta == 0))
})

test_that("initialize: init_theta = 'pl' sets theta to PL estimate", {
  s <- MRFPseudoBayes$new(z_small, R, fam, init_theta = "pl")
  pl_est <- mrf2d::smr_array(mrf2d::fit_pl(z_small, R, fam)$theta, fam)
  expect_equal(s$theta, pl_est)
})

test_that("initialize: numeric init_theta is accepted", {
  fdim <- length(mrf2d::smr_stat(z_small, R, fam))
  theta0 <- rep(0.5, fdim)
  s <- MRFPseudoBayes$new(z_small, R, fam, init_theta = theta0)
  expect_equal(s$theta, theta0)
})

test_that("initialize: invalid string init_theta throws error", {
  expect_error(
    MRFPseudoBayes$new(z_small, R, fam, init_theta = "bad"),
    "zero.*pl|pl.*zero"
  )
})

test_that("initialize: z must be a matrix", {
  expect_error(
    MRFPseudoBayes$new(as.vector(z_small), R, fam),
    "matrix"
  )
})

test_that("initialize: z with non-sequential labels throws error", {
  z_bad <- z_small
  z_bad[z_bad == 1] <- 5   # skips 1, 2, 3, 4
  expect_error(
    MRFPseudoBayes$new(z_bad, R, fam),
    "0, 1, \\.\\.\\., C"
  )
})

test_that("initialize: active bindings expose construction args", {
  s <- MRFPseudoBayes$new(z_small, R, fam, sdprior = 5, sdkernel = 0.1)
  expect_identical(s$z, z_small)
  expect_equal(s$family, fam)
  expect_equal(s$sdprior, 5)
  expect_equal(s$sdkernel, 0.1)
})

# ---------------------------------------------------------------------------
# Active bindings: sdprior / sdkernel are read-only
# ---------------------------------------------------------------------------

test_that("sdprior and sdkernel are read-only", {
  s <- MRFPseudoBayes$new(z_small, R, fam, sdprior = 5, sdkernel = 0.1)
  expect_equal(s$sdprior, 5)
  expect_equal(s$sdkernel, 0.1)
  expect_error(s$sdprior  <- 2)
  expect_error(s$sdkernel <- 0.2)
})

# ---------------------------------------------------------------------------
# Active binding: z is writable only with same dimensions
# ---------------------------------------------------------------------------

test_that("z can be replaced with a field of the same dimensions", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  z_new <- matrix(sample(0:2, prod(dim(z_small)), replace = TRUE),
                  nrow(z_small), ncol(z_small))
  expect_no_error(s$z <- z_new)
  expect_identical(s$z, z_new)
})

test_that("z cannot be replaced with different dimensions", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  z_wrong <- matrix(sample(0:2, 10 * 10, replace = TRUE), 10, 10)
  expect_error(s$z <- z_wrong, "Dimensions")
})

test_that("z cannot be replaced with a non-matrix", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  expect_error(s$z <- as.vector(z_small), "matrix")
})

# ---------------------------------------------------------------------------
# run()
# ---------------------------------------------------------------------------

test_that("run() returns self invisibly (chainable)", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  ret <- s$run(10, verbose = FALSE)
  expect_identical(ret, s)
})

test_that("run() increments n_samples", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(50, verbose = FALSE)
  expect_equal(s$n_samples, 50L)
})

test_that("run() called twice appends samples", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(30, verbose = FALSE)
  s$run(20, verbose = FALSE)
  expect_equal(s$n_samples, 50L)
  expect_equal(max(s$samples$t), 50L)
})

test_that("run() updates the chain state (theta)", {
  set.seed(1)
  s <- MRFPseudoBayes$new(z_small, R, fam, init_theta = "zero")
  theta_before <- s$theta
  s$run(200, verbose = FALSE)
  # With 200 steps the chain will almost certainly have moved
  expect_false(isTRUE(all.equal(s$theta, theta_before)))
})

# ---------------------------------------------------------------------------
# samples active binding
# ---------------------------------------------------------------------------

test_that("samples returns a tibble with required columns", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(20, verbose = FALSE)
  df <- s$samples
  expect_s3_class(df, "tbl_df")
  expect_named(df, c("t", "position", "interaction", "value"),
               ignore.order = FALSE)
})

test_that("samples has nrow == n_samples * fdim", {
  fdim <- length(mrf2d::smr_stat(z_small, R, fam))
  nrun <- 30L
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(nrun, verbose = FALSE)
  expect_equal(nrow(s$samples), nrun * fdim)
})

test_that("samples t column runs from 1 to n_samples", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(25, verbose = FALSE)
  expect_equal(sort(unique(s$samples$t)), seq_len(25L))
})

# ---------------------------------------------------------------------------
# summary()
# ---------------------------------------------------------------------------

test_that("summary() errors when no samples exist", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  expect_error(s$summary(), "No samples")
})

test_that("summary() returns a tibble with required columns", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(100, verbose = FALSE)
  stts <- s$summary()
  expect_s3_class(stts, "tbl_df")
  expect_true(all(c("position", "interaction", "q025", "mean", "q975", "sd")
                  %in% names(stts)))
})

test_that("summary() q025 <= mean <= q975", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(100, verbose = FALSE)
  stts <- s$summary()
  expect_true(all(stts$q025 <= stts$mean))
  expect_true(all(stts$mean <= stts$q975))
})

test_that("summary() burnin as fraction discards correct rows", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(100, verbose = FALSE)
  # Just check it doesn't error and returns same structure
  stts_frac <- s$summary(burnin = 0.5)
  stts_abs  <- s$summary(burnin = 50)
  expect_equal(stts_frac, stts_abs)
})

# ---------------------------------------------------------------------------
# plot()
# ---------------------------------------------------------------------------

test_that("plot() errors when no samples exist", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  expect_error(s$plot(), "No samples")
})

test_that("plot() returns a ggplot object", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(50, verbose = FALSE)
  p <- s$plot()
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# print()
# ---------------------------------------------------------------------------

test_that("print() runs without error", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  expect_no_error(s$print())
})

test_that("print() mentions family and n_samples", {
  s <- MRFPseudoBayes$new(z_small, R, fam)
  s$run(10, verbose = FALSE)
  out <- capture.output(s$print())
  expect_true(any(grepl(fam, out)))
  expect_true(any(grepl("10", out)))
})

# ---------------------------------------------------------------------------
# Multi-position / richer family smoke test
# ---------------------------------------------------------------------------

test_that("run() works with mrfi(2) and 'oneeach' family", {
  s <- MRFPseudoBayes$new(z_small, mrfi(2), "oneeach")
  s$run(20, verbose = FALSE)
  expect_equal(s$n_samples, 20L)
  expect_s3_class(s$samples, "tbl_df")
})
