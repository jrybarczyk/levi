### Tests
test_that("ok!", {
    base <- readExpColumn(a="Test-Control")
    expect_is(base, "list")
})

test_that("ok!", {
    expect_output(LEVIui(), "Parameter must be TRUE or FALSE")
})

