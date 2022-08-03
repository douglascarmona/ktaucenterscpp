test_that("rho_opt works", {
    t <- c(0.5591997, -0.1897360, 1.1834565, -5.9957109)
    expected_output <- c(0.048108, 0.005538, 0.215472, 1.000000)
    expect_equal(round(rho_opt(t = t, c = 1), 6), expected_output)
})