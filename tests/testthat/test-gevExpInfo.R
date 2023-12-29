# Check that the calculation of the expected information for the GEV
# distribution agrees with mev::gev.infomat()

sigma <- 1
mevMat <- cbind(
c(1/sigma^2, -0.422784335098467/sigma^2, 0.41184033042644/sigma),
c(-0.422784335098467/sigma^2, 1.82368066085288/sigma^2, 0.332484907160274/sigma),
c(0.41184033042644/sigma, 0.332484907160274/sigma, 2.42360605517703)
)

# Based on the constants
evilsMat1 <- cbind(c(imm(1, 0), ims0Constant, imx0Constant),
                   c(ims0Constant, iss0Constant, isx0Constant),
                   c(imx0Constant, isx0Constant, ixx0Constant))

# Based on calling the functions with xi = 0
evilsMat2 <- cbind(c(imm(1, 0), ims(1, 0), imx(1, 0)),
                   c(ims(1, 0), iss(1, 0), isx(1, 0)),
                   c(imx(1, 0), isx(1, 0), ixx(0)))

test_that("Expected Information constants agree with mev::gev.infomat()", {
  testthat::expect_equal(mevMat, evilsMat1)
})

test_that("Expected Information constants agree with mev::gev.infomat()", {
  testthat::expect_equal(mevMat, evilsMat2)
})
