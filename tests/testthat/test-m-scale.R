test_that("normal_consistency_constants works", {
    expect_equal(normal_consistency_constants(1), 0.404629)
    expect_equal(normal_consistency_constants(10), 1.739537)
    expect_equal(normal_consistency_constants(400), 11.09393)
})