library(mrf2d)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------
set.seed(42)
nr <- 15L; nc <- 15L
z_small <- matrix(sample(0:1, nr * nc, replace = TRUE), nr, nc)  # binary C=1
R_max   <- mrfi(1)          # 2 positions: (1,0) and (0,1)
fam     <- "oneeach"        # dpg=1 for C=1 → fdim=2
npos    <- length(R_max)    # 2

new_rj <- function(...) {
  MRFRJPseudoBayes$new(z_small, R_max, fam, ...)
}

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

test_that("initialize returns an R6 object of class MRFRJPseudoBayes", {
  s <- new_rj()
  expect_s3_class(s, "MRFRJPseudoBayes")
  expect_true(R6::is.R6(s))
})

test_that("MRFRJPseudoBayes inherits from MRFPseudoBayes", {
  s <- new_rj()
  expect_s3_class(s, "MRFPseudoBayes")
})

test_that("initialize: n_samples starts at 0 and samples/samples_mrfi are NULL", {
  s <- new_rj()
  expect_equal(s$n_samples, 0L)
  expect_null(s$samples)
  expect_null(s$samples_mrfi)
})

test_that("initialize: init_included = 'zero' sets all positions inactive", {
  s <- new_rj(init_included = "zero")
  expect_true(all(!s$included))
  expect_equal(s$n_included, 0L)
})

test_that("initialize: init_included = 'full' sets all positions active", {
  s <- new_rj(init_included = "full")
  expect_true(all(s$included))
  expect_equal(s$n_included, npos)
})

test_that("initialize: logical vector init_included is accepted", {
  inc <- c(TRUE, FALSE)
  s   <- new_rj(init_included = inc)
  expect_equal(s$included, inc)
  expect_equal(s$n_included, 1L)
})

test_that("initialize: invalid init_included string throws error", {
  expect_error(new_rj(init_included = "bad"), "zero.*full|full.*zero")
})

test_that("initialize: wrong-length init_included throws error", {
  expect_error(new_rj(init_included = c(TRUE, FALSE, FALSE)))  # length 3 != npos=2
})

test_that("initialize: active bindings expose scalar parameters", {
  s <- new_rj(sdprior = 2, sdkernel = 0.01, sdbirth = 0.02,
              logpenalty = 5)
  expect_equal(s$sdprior,    2)
  expect_equal(s$sdkernel,   0.01)
  expect_equal(s$sdbirth,    0.02)
  expect_equal(s$logpenalty, 5)
})

test_that("initialize: default logpenalty equals log(n*m)", {
  s <- new_rj()
  expect_equal(s$logpenalty, log(nr * nc))
})

test_that("initialize: theta zeros match init_included='zero'", {
  s <- new_rj(init_included = "zero")
  expect_true(all(s$theta == 0))
})

test_that("initialize: init_theta='pl' restricts to active positions", {
  inc <- c(TRUE, FALSE)
  s   <- new_rj(init_included = inc, init_theta = "pl")
  # Inactive position (index 2) must have theta = 0
  dpg <- length(s$theta) / npos   # dpg = 1
  inactive_idx <- ((2L - 1L) * dpg + 1L):(2L * dpg)  # index 2
  expect_true(all(s$theta[inactive_idx] == 0))
})

# ---------------------------------------------------------------------------
# included is read-only
# ---------------------------------------------------------------------------

test_that("included active binding is read-only", {
  s <- new_rj()
  expect_error(s$included <- rep(TRUE, npos))
})

# ---------------------------------------------------------------------------
# run()
# ---------------------------------------------------------------------------

test_that("run() returns self invisibly (chainable)", {
  s   <- new_rj()
  ret <- s$run(10, verbose = FALSE)
  expect_identical(ret, s)
})

test_that("run() increments n_samples", {
  s <- new_rj()
  s$run(30, verbose = FALSE)
  expect_equal(s$n_samples, 30L)
})

test_that("run() appends on repeated calls", {
  s <- new_rj()
  s$run(20, verbose = FALSE)
  s$run(10, verbose = FALSE)
  expect_equal(s$n_samples, 30L)
})

test_that("run() starting from 'full' can reduce active set", {
  set.seed(7)
  # Large logpenalty strongly discourages large models; birth/death both non-zero
  s <- new_rj(init_included = "full",
              logpenalty   = 15,
              kernel_probs = c(1, 1, 5, 5, 5))
  s$run(500, verbose = FALSE)
  # With high penalty and random z (no structure), smaller model is preferred
  expect_lt(s$n_included, npos)
})

test_that("run() starting from 'zero' can grow active set", {
  set.seed(13)
  # Negative logpenalty rewards larger models, so added positions stay active
  s <- new_rj(init_included = "zero",
              logpenalty   = -10,
              kernel_probs = c(1, 1, 3, 3, 3))
  s$run(200, verbose = FALSE)
  expect_gt(s$n_included, 0L)
})

# ---------------------------------------------------------------------------
# samples active binding (filtered zeros)
# ---------------------------------------------------------------------------

test_that("samples returns NULL before run()", {
  expect_null(new_rj()$samples)
})

test_that("samples returns tibble with required columns after run()", {
  s <- new_rj(init_included = "full")
  s$run(20, verbose = FALSE)
  df <- s$samples
  expect_s3_class(df, "tbl_df")
  expect_named(df, c("t", "position", "interaction", "value"),
               ignore.order = FALSE)
})

test_that("samples contains only nonzero values", {
  s <- new_rj(init_included = "full")
  s$run(20, verbose = FALSE)
  expect_true(all(s$samples$value != 0))
})

test_that("samples has t column bounded by n_samples", {
  s <- new_rj(init_included = "full")
  s$run(25, verbose = FALSE)
  expect_true(all(s$samples$t >= 1L & s$samples$t <= 25L))
})

# ---------------------------------------------------------------------------
# samples_mrfi active binding
# ---------------------------------------------------------------------------

test_that("samples_mrfi returns NULL before run()", {
  expect_null(new_rj()$samples_mrfi)
})

test_that("samples_mrfi returns tibble with required columns", {
  s <- new_rj()
  s$run(15, verbose = FALSE)
  df <- s$samples_mrfi
  expect_s3_class(df, "tbl_df")
  expect_named(df, c("t", "position", "value"), ignore.order = FALSE)
})

test_that("samples_mrfi has npos rows per iteration", {
  nrun <- 10L
  s    <- new_rj()
  s$run(nrun, verbose = FALSE)
  expect_equal(nrow(s$samples_mrfi), nrun * npos)
})

test_that("samples_mrfi value column is logical", {
  s <- new_rj()
  s$run(10, verbose = FALSE)
  expect_type(s$samples_mrfi$value, "logical")
})

# ---------------------------------------------------------------------------
# summary()
# ---------------------------------------------------------------------------

test_that("summary() errors when no samples exist", {
  s <- new_rj()
  expect_error(s$summary(), "No samples")
})

test_that("summary() returns tibble with required columns", {
  s <- new_rj(init_included = "full")
  s$run(100, verbose = FALSE)
  st <- s$summary()
  expect_s3_class(st, "tbl_df")
  expect_true(all(c("position", "interaction", "prob",
                    "q025", "mean", "q975", "sd") %in% names(st)))
})

test_that("summary() inclusion probabilities are in [0, 1]", {
  s <- new_rj(init_included = "full")
  s$run(100, verbose = FALSE)
  st <- s$summary()
  expect_true(all(st$prob >= 0 & st$prob <= 1))
})

test_that("summary() contains all positions including never-active ones", {
  s <- new_rj(init_included = "zero")  # no position ever active
  s$run(20, verbose = FALSE)
  st <- s$summary()
  expect_equal(nrow(st), private_fdim <- length(s$theta))  # fdim rows total
  expect_true(all(st$prob == 0))
  expect_true(all(is.na(st$q025)))
})

test_that("summary() credible bounds are ordered for active positions", {
  # Negative logpenalty ensures positions stay active throughout the run
  s <- new_rj(init_included = "full", logpenalty = -10)
  s$run(100, verbose = FALSE)
  st <- s$summary()
  active <- st[!is.na(st$q025), ]
  expect_gt(nrow(active), 0L)  # at least one active position expected
  expect_true(all(active$q025 <= active$mean))
  expect_true(all(active$mean <= active$q975))
})

# ---------------------------------------------------------------------------
# plot()
# ---------------------------------------------------------------------------

test_that("plot() returns a ggplot", {
  s <- new_rj(init_included = "full")
  s$run(100, verbose = FALSE)
  p <- s$plot()
  expect_s3_class(p, "gg")
})

test_that("plot() errors when no samples exist", {
  s <- new_rj()
  expect_error(s$plot(), "No samples")
})

# ---------------------------------------------------------------------------
# plot_mrfi()
# ---------------------------------------------------------------------------

test_that("plot_mrfi() returns a ggplot after run()", {
  s <- new_rj()
  s$run(30, verbose = FALSE)
  p <- s$plot_mrfi()
  expect_s3_class(p, "gg")
})

test_that("plot_mrfi() errors when no samples exist", {
  s <- new_rj()
  expect_error(s$plot_mrfi(), "No samples")
})

# ---------------------------------------------------------------------------
# print()
# ---------------------------------------------------------------------------

test_that("print() mentions class name", {
  s   <- new_rj()
  out <- capture.output(s$print())
  expect_true(any(grepl("MRFRJPseudoBayes", out)))
})

test_that("print() shows active/total positions", {
  s   <- new_rj(init_included = c(TRUE, FALSE))
  out <- capture.output(s$print())
  expect_true(any(grepl("1/2", out)))
})

# ---------------------------------------------------------------------------
# Inherited methods still work
# ---------------------------------------------------------------------------

test_that("inherited z, mrfi, family, sdprior, sdkernel bindings are accessible", {
  s <- new_rj(sdprior = 3, sdkernel = 0.02)
  expect_identical(s$z, z_small)
  expect_equal(s$family, fam)
  expect_equal(s$sdprior, 3)
  expect_equal(s$sdkernel, 0.02)
})
