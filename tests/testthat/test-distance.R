test_that("distance_to_centers works", {
    set.seed(10)
    Z <- rnorm(400)
    mues <- rep(c(0.0, 4.0), 200)
    data <-  matrix(Z + mues, ncol = 2)
    row_odd <- seq_len(nrow(data)) %% 2
    centers <- rbind(c(0.0, 0.0), c(4.0, 4.0))

    dists <- distance_to_centers(data, centers)
    distances_min <- dists$min_distance
    cluster <- dists$membership

    expect_equal(cluster[row_odd == 1], rep(1, 100))
    expect_equal(cluster[row_odd == 0], rep(2, 100))
    expect_equal(round(distances_min[1], 6), 1.215658)
    expect_equal(round(distances_min[200], 6), 1.660615)
})