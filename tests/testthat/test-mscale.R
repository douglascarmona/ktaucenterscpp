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
