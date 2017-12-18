context("Funct")

test_that("Linear Solver works",{
  A1 <- c(2,3,4,9)
  b1 <- c(6,15)
  A2 <- c(1,3,-2,3,5,6,2,4,3)
  b2 <- c(5,7,8)
  A3 <- c(1,1,-2,1,3,-1,2,-1,1,2,2,-3,1,3,-3,-1,2,1,5,2,-1,-1,2,1,-3,-1,2,3,1,3,4,3,1,-6,-3,-2)
  b3 <- c(4,20,-15,-3,16,-27)
  Exp_value1 <- matrix(c(1.5,1),nrow = 2,ncol=1,byrow = TRUE)
  Exp_value2 <- matrix(c(-15,8,2),nrow = 3,ncol=1,byrow = TRUE)
  Exp_value3 <- matrix(c(1/3,-430/99,313/99,104/99,142/33,-37/99),nrow = 6,ncol=1,byrow = TRUE)
  expect_equal(Funct(A1,b1,2),Exp_value1)
  expect_equal(Funct(A2,b2,2),Exp_value2)       #In case the values don't match then we get an error along with the difference in the value of the vector
  expect_equal(Funct(A3,b3,2),Exp_value3)
})

