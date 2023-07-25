test_that("mrfrj works", {
  z <- mrf2d::Z_potts
  ll <- llapprox(z, mrfi(2), family = "oneeach", method = "pseudo")
  result <- mrfrj(z = z,llapprox = ll)
  expect_s3_class(result, "mrfbayes_out", TRUE)
})
