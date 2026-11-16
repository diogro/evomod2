test_that("autonomy is 1 along eigenvectors and matches the definition", {
  G <- matrix(c(2, 0.8, 0.3, 0.8, 1, 0.1, 0.3, 0.1, 0.5), 3)
  v <- eigen(G)$vectors
  expect_equal(autonomy(G, v), rep(1, 3))
  b <- c(1, 2, 0)
  bu <- b / sqrt(sum(b^2))
  expect_equal(autonomy(G, b), 1 / (drop(t(bu) %*% solve(G) %*% bu) * drop(t(bu) %*% G %*% bu)))
  expect_length(autonomy(G, cbind(b, 1:3)), 2)
  set.seed(1)
  a <- mean_autonomy(G, 500)
  expect_true(a > 0 && a < 1)
  expect_equal(sum(eigen_variance(G)), 1)
  expect_true(all(diff(eigen_variance(G)) <= 0))
})

test_that("autonomy agrees with evolqg", {
  skip_if_not_installed("evolqg")
  G <- matrix(c(2, 0.8, 0.3, 0.8, 1, 0.1, 0.3, 0.1, 0.5), 3)
  b <- cbind(c(1, 2, 0), c(0, 1, 1))
  b <- sweep(b, 2, sqrt(colSums(b^2)), "/")
  expect_equal(unname(autonomy(G, b)), unname(evolqg::Autonomy(G, b)))
})
