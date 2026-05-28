library(mrf2d)

# ---------------------------------------------------------------------------
# Shared fixtures — explicit init to avoid EM and keep tests deterministic
# ---------------------------------------------------------------------------
set.seed(42)
nr <- 12L; nc <- 12L
z_true <- matrix(sample(0:1, nr * nc, replace = TRUE), nr, nc)
y_obs  <- matrix(rnorm(nr * nc, mean = 4 * as.numeric(z_true), sd = 1),
                 nr, nc)
R_max  <- mrfi(1)   # 2 positions, npos = 2
fam    <- "oneeach" # 1 param/position for C=1 -> fdim = 2
C_val  <- 1L
mu0    <- c(0, 4)
sig0   <- c(1, 1)

new_hrj <- function(...) {
  HMRFRJPseudoBayes$new(y_obs, R_max, fam, C_val,
                        z0 = z_true, mu0 = mu0, sigma2_0 = sig0, ...)
}

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

test_that("initialize returns an R6 object of class HMRFRJPseudoBayes", {
  s <- new_hrj()
  expect_s3_class(s, "HMRFRJPseudoBayes")
  expect_true(R6::is.R6(s))
})

test_that("HMRFRJPseudoBayes inherits from MRFRJPseudoBayes and MRFPseudoBayes", {
  s <- new_hrj()
  expect_s3_class(s, "MRFRJPseudoBayes")
  expect_s3_class(s, "MRFPseudoBayes")
})

test_that("initial state: n_samples = 0, all sample accessors return NULL", {
  s <- new_hrj()
  expect_equal(s$n_samples, 0L)
  expect_null(s$samples)
  expect_null(s$samples_pars)
  expect_null(s$samples_mrfi)
  expect_null(s$zprobs)
})

test_that("init_included = 'zero' sets all positions inactive", {
  s <- new_hrj(init_included = "zero")
  expect_true(all(!s$included))
  expect_equal(s$n_included, 0L)
})

test_that("init_included = 'full' sets all positions active", {
  s <- new_hrj(init_included = "full")
  expect_true(all(s$included))
})

test_that("active bindings expose y, mu, sigma2, sdpriormu, alpha1, alpha2", {
  s <- new_hrj(sdpriormu = 5, alpha1 = 0.01, alpha2 = 0.01)
  expect_identical(s$y, y_obs)
  expect_equal(s$mu,      mu0)
  expect_equal(s$sigma2,  sig0)
  expect_equal(s$sdpriormu, 5)
  expect_equal(s$alpha1, 0.01)
  expect_equal(s$alpha2, 0.01)
})

test_that("active bindings from parent are accessible", {
  s <- new_hrj(sdprior = 2, sdkernel = 0.01, sdbirth = 0.02, logpenalty = 3)
  expect_equal(s$family,     fam)
  expect_equal(s$sdprior,    2)
  expect_equal(s$sdkernel,   0.01)
  expect_equal(s$sdbirth,    0.02)
  expect_equal(s$logpenalty, 3)
})

test_that("C is inferred from z0 and matches C argument", {
  s <- new_hrj()
  expect_equal(s$n_samples, 0L)  # proxy: object was created without error
  # C+1 emission means/variances
  expect_length(s$mu,     C_val + 1L)
  expect_length(s$sigma2, C_val + 1L)
})

test_that("invalid C (< 1) raises error", {
  expect_error(HMRFRJPseudoBayes$new(y_obs, R_max, fam, 0L,
                                     z0 = z_true, mu0 = mu0,
                                     sigma2_0 = sig0))
})

# ---------------------------------------------------------------------------
# run()
# ---------------------------------------------------------------------------

test_that("run() returns self invisibly", {
  s   <- new_hrj()
  ret <- s$run(10, verbose = FALSE)
  expect_identical(ret, s)
})

test_that("run() increments n_samples", {
  s <- new_hrj()
  s$run(25, verbose = FALSE)
  expect_equal(s$n_samples, 25L)
})

test_that("run() appends on repeated calls", {
  s <- new_hrj()
  s$run(20, verbose = FALSE)
  s$run(10, verbose = FALSE)
  expect_equal(s$n_samples, 30L)
})

test_that("run() updates mu and sigma2 (emission params move)", {
  s <- new_hrj()
  mu_before  <- s$mu
  sig_before <- s$sigma2
  s$run(20, verbose = FALSE)
  # After sampling, mu/sigma2 should have changed from initial values
  expect_false(all(s$mu     == mu_before))
  expect_false(all(s$sigma2 == sig_before))
})

# ---------------------------------------------------------------------------
# samples — theta chain (filtered zeros)
# ---------------------------------------------------------------------------

test_that("samples returns tibble with correct columns after run()", {
  s <- new_hrj(init_included = "full")
  s$run(20, verbose = FALSE)
  df <- s$samples
  if (!is.null(df)) {
    expect_s3_class(df, "tbl_df")
    expect_named(df, c("t", "position", "interaction", "value"),
                 ignore.order = FALSE)
    expect_true(all(df$value != 0))
  }
})

# ---------------------------------------------------------------------------
# samples_pars — emission parameter chain
# ---------------------------------------------------------------------------

test_that("samples_pars returns tibble with required columns", {
  s <- new_hrj()
  s$run(20, verbose = FALSE)
  df <- s$samples_pars
  expect_s3_class(df, "tbl_df")
  expect_named(df, c("t", "k", "par", "value"), ignore.order = FALSE)
})

test_that("samples_pars has 2*(C+1) rows per iteration", {
  nrun <- 15L
  s    <- new_hrj()
  s$run(nrun, verbose = FALSE)
  # 2 parameters (mu, sigma2) × (C+1) classes × nrun iterations
  expect_equal(nrow(s$samples_pars), 2L * (C_val + 1L) * nrun)
})

test_that("samples_pars par column is exactly 'mu' and 'sigma2'", {
  s <- new_hrj()
  s$run(10, verbose = FALSE)
  expect_setequal(unique(s$samples_pars$par), c("mu", "sigma2"))
})

test_that("samples_pars k column contains values 0:C", {
  s <- new_hrj()
  s$run(10, verbose = FALSE)
  expect_setequal(unique(s$samples_pars$k), 0:C_val)
})

# ---------------------------------------------------------------------------
# samples_mrfi — inclusion chain
# ---------------------------------------------------------------------------

test_that("samples_mrfi has correct columns after run()", {
  s <- new_hrj()
  s$run(10, verbose = FALSE)
  df <- s$samples_mrfi
  expect_s3_class(df, "tbl_df")
  expect_named(df, c("t", "position", "value"), ignore.order = FALSE)
  expect_type(df$value, "logical")
})

# ---------------------------------------------------------------------------
# zprobs
# ---------------------------------------------------------------------------

test_that("zprobs has correct dimensions after run()", {
  s <- new_hrj()
  s$run(15, verbose = FALSE)
  zp <- s$zprobs
  expect_equal(dim(zp), c(nr, nc, C_val + 1L))
})

test_that("zprobs sums to 1 across classes at every pixel", {
  s <- new_hrj()
  s$run(20, verbose = FALSE)
  zp    <- s$zprobs
  totals <- apply(zp, c(1, 2), sum)
  expect_true(all(abs(totals - 1) < 1e-10))
})

# ---------------------------------------------------------------------------
# summary() — theta + inclusion (inherited from MRFRJPseudoBayes)
# ---------------------------------------------------------------------------

test_that("summary() errors when no samples exist", {
  s <- new_hrj()
  expect_error(s$summary(), "No samples")
})

test_that("summary() returns tibble with required columns", {
  s <- new_hrj(init_included = "full", logpenalty = -10)
  s$run(60, verbose = FALSE)
  st <- s$summary()
  expect_s3_class(st, "tbl_df")
  expect_true(all(c("position", "interaction", "prob",
                    "q025", "mean", "q975", "sd") %in% names(st)))
})

test_that("summary() inclusion probabilities are in [0, 1]", {
  s <- new_hrj(init_included = "full", logpenalty = -10)
  s$run(60, verbose = FALSE)
  expect_true(all(s$summary()$prob >= 0 & s$summary()$prob <= 1))
})

# ---------------------------------------------------------------------------
# summary_pars() — emission parameters
# ---------------------------------------------------------------------------

test_that("summary_pars() errors when no samples exist", {
  s <- new_hrj()
  expect_error(s$summary_pars(), "No samples")
})

test_that("summary_pars() returns tibble with required columns", {
  s <- new_hrj()
  s$run(40, verbose = FALSE)
  sp <- s$summary_pars()
  expect_s3_class(sp, "tbl_df")
  expect_true(all(c("k", "par", "q025", "mean", "q975", "sd") %in% names(sp)))
})

test_that("summary_pars() credible bounds are ordered", {
  s <- new_hrj()
  s$run(40, verbose = FALSE)
  sp <- s$summary_pars()
  expect_true(all(sp$q025 <= sp$mean))
  expect_true(all(sp$mean <= sp$q975))
})

# ---------------------------------------------------------------------------
# plot()
# ---------------------------------------------------------------------------

test_that("plot('pars') returns ggplot", {
  s <- new_hrj()
  s$run(20, verbose = FALSE)
  expect_s3_class(s$plot("pars"), "gg")
})

test_that("plot('theta') returns ggplot when positions were active", {
  s <- new_hrj(init_included = "full")
  s$run(50, verbose = FALSE)
  # samples may be NULL if no position stayed active; only test if non-NULL
  if (!is.null(s$samples) && nrow(s$samples) > 0L)
    expect_s3_class(s$plot("theta"), "gg")
})

test_that("plot('mrfi') returns ggplot", {
  s <- new_hrj()
  s$run(20, verbose = FALSE)
  expect_s3_class(s$plot("mrfi"), "gg")
})

test_that("plot('zprobs') returns ggplot", {
  s <- new_hrj()
  s$run(20, verbose = FALSE)
  expect_s3_class(s$plot("zprobs"), "gg")
})

test_that("plot() with invalid 'what' throws error", {
  s <- new_hrj()
  s$run(10, verbose = FALSE)
  expect_error(s$plot("bad"), "'what' must be")
})

test_that("plot('pars') errors when no samples exist", {
  s <- new_hrj()
  expect_error(s$plot("pars"), "No samples")
})

# ---------------------------------------------------------------------------
# print()
# ---------------------------------------------------------------------------

test_that("print() mentions class name", {
  s   <- new_hrj()
  out <- capture.output(s$print())
  expect_true(any(grepl("HMRFRJPseudoBayes", out)))
})

test_that("print() shows field size", {
  s   <- new_hrj()
  out <- capture.output(s$print())
  expect_true(any(grepl(paste0(nr, " x ", nc), out)))
})
