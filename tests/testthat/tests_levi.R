### Tests
test_that("ok!", {
    base <- readExpColumn(a="Test-Control")
    expect_is(base, "list")
})

test_that("ok!", {
    expect_output(LEVIui(), "Parameter must be TRUE or FALSE")
})


# Teste a função matrix_entrada
test_that("matrix_entrada works correctly", {
  coordPiso <- matrix(c(1, 0, 1, 0, 1, 1, 0, 1, 0), nrow = 3, ncol = 3)
  SignalOut <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), nrow = 9, ncol = 1)
  signalExp <- matrix(c(0.5, 0.4, 0.3, 0.2, 0.1, 0.6, 0.7, 0.8, 0.9), nrow = 9, ncol = 1)
  signalCtrl <- matrix(c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1), nrow = 9, ncol = 1)
  coord <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), nrow = 9, ncol = 2)
  resolutionValue <- 3
  increase <- 0.1
  zoomValue <- 0.0
  numberCoord <- 9
  
  result <- matrix_entrada(coordPiso, SignalOut, signalExp, signalCtrl, coord, resolutionValue, increase, zoomValue, numberCoord)
  
  expect_equal(dim(result$m1), c(18, 5)) # numberCoord + (resolutionValue * resolutionValue)
  expect_equal(result$m3, 5) # número esperado de entradas na matriz
  
  # Verifique se as primeiras entradas correspondem às coordenadas fornecidas
  expect_equal(result$m1[1:9, 1], c(0.0,0.0,0.1,0.1,0.2,0.0,0.0,0.0,0.0))
  expect_equal(result$m1[1:9, 2], c(0.1,0.2,0.0,0.2,0.2,0.0,0.0,0.0,0.0))
  expect_equal(result$m1[1:9, 3], c(0.0,0.0,0.0,0.1,0.0,0.0,0.0,0.0,0.0))
  expect_equal(result$m1[1:9, 4], c(0.0,0.0,0.0,0.5,0.0,0.0,0.0,0.0,0.0))
  expect_equal(result$m1[1:9, 5], c(0.0,0.0,0.0,0.9,0.0,0.0,0.0,0.0,0.0))
  
  # Verifique se as últimas entradas são zeros (com base em coordPiso)
  expect_equal(result$m1[10:12, 3:5], matrix(0, ncol = 3, nrow = 3))
})