library(mrf2d)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
# Small synthetic HMRF dataset (C = 1, two Gaussian classes)
set.seed(7)
nr <- 12L; nc <- 12L
z_true <- matrix(sample(0:1, nr * nc, replace = TRUE), nr, nc)
y_obs  <- matrix(
  ifelse(as.vector(z_true) == 0L,
         rnorm(nr * nc, mean = 0, sd = 1),
         rnorm(nr * nc, mean = 4, sd = 1)),
  nr, nc
)
R   <- mrfi(1)
fam <- "onepar"
C   <- 1L

# Known-good initial values (bypass EM in all tests)
mu0      <- c(0, 4)
sigma2_0 <- c(1, 1)

new_sampler <- function(...) {
  HMRFPseudoBayes$new(
    y = y_obs, mrfi = R, family = fam, C = C,
    z0 = z_true, mu0 = mu0, sigma2_0 = sigma2_0,
    verbose = FALSE, ...
  )
}

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

test_that("initialize returns an R6 object of class HMRFPseudoBayes", {
  s <- new_sampler()
  expect_s3_class(s, "HMRFPseudoBayes")
  expect_true(R6::is.R6(s))
})

test_that("HMRFPseudoBayes inherits from MRFPseudoBayes", {
  s <- new_sampler()
  expect_s3_class(s, "MRFPseudoBayes")
})

test_that("initialize: n_samples starts at 0 and samples/samples_pars are NULL", {
  s <- new_sampler()
  expect_equal(s$n_samples, 0L)
  expect_null(s$samples)
  expect_null(s$samples_pars)
})

test_that("initialize: active bindings expose construction arguments", {
  s <- new_sampler(sdprior = 5, sdkernel = 0.01,
                   sdpriormu = 8, alpha1 = 0.1, alpha2 = 0.2)
  expect_identical(s$y, y_obs)
  expect_equal(s$mu, mu0)
  expect_equal(s$sigma2, sigma2_0)
  expect_equal(s$sdpriormu, 8)
  expect_equal(s$alpha1, 0.1)
  expect_equal(s$alpha2, 0.2)
  expect_equal(s$sdprior, 5)
  expect_equal(s$sdkernel, 0.01)
})

test_that("initialize: y must be a numeric matrix", {
  expect_error(
    HMRFPseudoBayes$new(as.vector(y_obs), R, fam, C,
                        z0 = z_true, mu0 = mu0, sigma2_0 = sigma2_0),
    "matrix"
  )
})

test_that("initialize: C < 1 throws error", {
  expect_error(
    HMRFPseudoBayes$new(y_obs, R, fam, C = 0L,
                        z0 = z_true, mu0 = mu0, sigma2_0 = sigma2_0,
                        verbose = FALSE)
  )
})

test_that("initialize: sigma2_0 must be positive", {
  expect_error(
    HMRFPseudoBayes$new(y_obs, R, fam, C,
                        z0 = z_true, mu0 = mu0, sigma2_0 = c(-1, 1),
                        verbose = FALSE)
  )
})

# ---------------------------------------------------------------------------
# run()
# ---------------------------------------------------------------------------

test_that("run() returns self invisibly (chainable)", {
  s   <- new_sampler()
  ret <- s$run(5, verbose = FALSE)
  expect_identical(ret, s)
})

test_that("run() increments n_samples for theta chain", {
  s <- new_sampler()
  s$run(20, verbose = FALSE)
  expect_equal(s$n_samples, 20L)
})

test_that("run() appends samples on repeated calls", {
  s <- new_sampler()
  s$run(15, verbose = FALSE)
  s$run(10, verbose = FALSE)
  expect_equal(s$n_samples, 25L)
  expect_equal(max(s$samples$t), 25L)
})

test_that("run() mu and sigma2 change after sampling", {
  set.seed(99)
  s      <- new_sampler()
  mu_b   <- s$mu
  sig_b  <- s$sigma2
  s$run(30, verbose = FALSE)
  expect_false(isTRUE(all.equal(s$mu,     mu_b)))
  expect_false(isTRUE(all.equal(s$sigma2, sig_b)))
})

test_that("run() latent field z changes after sampling", {
  set.seed(11)
  s    <- new_sampler()
  z_b  <- s$z
  s$run(10, verbose = FALSE)
  # z is categorical and small; very unlikely all pixels are identical
  expect_false(identical(s$z, z_b))
})

# ---------------------------------------------------------------------------
# samples (theta) — inherited active binding
# ---------------------------------------------------------------------------

test_that("samples returns tibble with required columns", {
  s <- new_sampler()
  s$run(15, verbose = FALSE)
  df <- s$samples
  expect_s3_class(df, "tbl_df")
  expect_named(df, c("t", "position", "interaction", "value"),
               ignore.order = FALSE)
})

# ---------------------------------------------------------------------------
# samples_pars active binding
# ---------------------------------------------------------------------------

test_that("samples_pars returns tibble with required columns", {
  s <- new_sampler()
  s$run(15, verbose = FALSE)
  df <- s$samples_pars
  expect_s3_class(df, "tbl_df")
  expect_named(df, c("t", "k", "par", "value"), ignore.order = FALSE)
})

test_that("samples_pars contains exactly 'mu' and 'sigma2' pars", {
  s <- new_sampler()
  s$run(10, verbose = FALSE)
  expect_setequal(unique(s$samples_pars$par), c("mu", "sigma2"))
})

test_that("samples_pars has correct number of rows", {
  nrun <- 20L
  s    <- new_sampler()
  s$run(nrun, verbose = FALSE)
  # nrun samples * (C+1) labels * 2 parameters (mu, sigma2)
  expect_equal(nrow(s$samples_pars), nrun * (C + 1L) * 2L)
})

test_that("samples_pars k values are 0 and 1 for C=1", {
  s <- new_sampler()
  s$run(10, verbose = FALSE)
  expect_setequal(unique(s$samples_pars$k), 0:1)
})

# ---------------------------------------------------------------------------
# zprobs active binding
# ---------------------------------------------------------------------------

test_that("zprobs is NULL before run()", {
  s <- new_sampler()
  expect_null(s$zprobs)
})

test_that("zprobs has correct dimensions", {
  s <- new_sampler()
  s$run(10, verbose = FALSE)
  zp <- s$zprobs
  expect_equal(dim(zp), c(nr, nc, C + 1L))
})

test_that("zprobs sums to 1 across labels at every pixel", {
  s <- new_sampler()
  s$run(20, verbose = FALSE)
  label_sums <- apply(s$zprobs, c(1, 2), sum)
  expect_true(all(abs(label_sums - 1) < 1e-10))
})

# ---------------------------------------------------------------------------
# summary_pars()
# ---------------------------------------------------------------------------

test_that("summary_pars() errors when no samples exist", {
  s <- new_sampler()
  expect_error(s$summary_pars(), "No samples")
})

test_that("summary_pars() returns tibble with required columns", {
  s <- new_sampler()
  s$run(60, verbose = FALSE)
  st <- s$summary_pars()
  expect_s3_class(st, "tbl_df")
  expect_true(all(c("k", "par", "q025", "mean", "q975", "sd") %in% names(st)))
})

test_that("summary_pars() credible bounds are ordered", {
  s <- new_sampler()
  s$run(60, verbose = FALSE)
  st <- s$summary_pars()
  expect_true(all(st$q025 <= st$mean))
  expect_true(all(st$mean <= st$q975))
})

# ---------------------------------------------------------------------------
# Inherited summary() still works for theta
# ---------------------------------------------------------------------------

test_that("summary() (theta) still works via inheritance", {
  s <- new_sampler()
  s$run(60, verbose = FALSE)
  st <- s$summary()
  expect_s3_class(st, "tbl_df")
  expect_true(all(c("position", "interaction", "q025", "mean", "q975", "sd")
                  %in% names(st)))
})

# ---------------------------------------------------------------------------
# plot() dispatch
# ---------------------------------------------------------------------------

test_that("plot(what='pars') returns a ggplot", {
  s <- new_sampler()
  s$run(10, verbose = FALSE)
  p <- s$plot(what = "pars")
  expect_s3_class(p, "gg")
})

test_that("plot(what='theta') returns a ggplot (delegates to parent)", {
  s <- new_sampler()
  s$run(10, verbose = FALSE)
  p <- s$plot(what = "theta")
  expect_s3_class(p, "gg")
})

test_that("plot(what='zprobs') returns a ggplot", {
  s <- new_sampler()
  s$run(10, verbose = FALSE)
  p <- s$plot(what = "zprobs")
  expect_s3_class(p, "gg")
})

test_that("plot() with invalid what throws error", {
  s <- new_sampler()
  s$run(10, verbose = FALSE)
  expect_error(s$plot(what = "bad"), "pars.*theta.*zprobs|theta.*pars.*zprobs")
})

# ---------------------------------------------------------------------------
# print()
# ---------------------------------------------------------------------------

test_that("print() does not error and mentions class name", {
  s   <- new_sampler()
  out <- capture.output(s$print())
  expect_true(any(grepl("HMRFPseudoBayes", out)))
})
