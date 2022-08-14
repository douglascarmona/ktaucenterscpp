test_that("normal_consistency_constants works", {
    expect_equal(normal_consistency_constants(1), 0.404629)
    expect_equal(normal_consistency_constants(10), 1.739537)
    expect_equal(normal_consistency_constants(400), 11.09393)
})

test_that("c1 works", {
    expect_equal(c1(), 1.0)
})

test_that("c2 works", {
    expect_equal(round(c2(1), 6), 2.9987)
    expect_equal(round(c2(10), 6), 1.028568)
})

test_that("mscale works", {
    u1 <- rep(0, 5)
    u2 <- c(
        -5.9957109, 3.7043719, 8.3375155, -4.3120109, -7.9069974,
        4.0211492, 0.5591997, 6.1587040, 9.1300025, -7.7909396
        )
    u3 <- c(
        -6.4290374, -9.7518368, -7.0511081, -3.4985002, -3.9629706,
        -2.4026255, -6.6018854, 2.3083618, 0.8768204, 0.5212299
        )
    u4 <- c(
        -7.136806, -12.296880, -16.190480, -14.629739, -11.182638,
        -11.644894, -10.332334, -5.966643, -15.903391, 1.307210
    )
    expect_equal(mscale(u = u1, c = 1.0, b = 0.5), 0)
    expect_equal(round(mscale(u = u2, c = 1.0, b = 0.5), 6), 3.444974)
    expect_equal(round(mscale(u = u3, c = 1.0, b = 0.5), 6), 2.417396)
    expect_equal(round(mscale(u = u4, c = 1.0, b = 0.5), 6), 6.316444)
})
