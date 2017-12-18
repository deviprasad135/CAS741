library(optR)
set.seed(2001)
#Gauss-Elimination
#Define Functions
########################################################
#Define a function to convert the coeffcient vector into a matrix form
Matrix_Conv_A <- function(A){
  n <- (length(A))^0.5
  return(matrix(A,nrow = n,ncol=n,byrow = TRUE))
}
#Define a function to convert the constant vector into a matrix form
Matrix_Conv_B <- function(b){
  n <- length(b)
  return(matrix(b,nrow=n,ncol=1,byrow = TRUE))
}
#Create a upper triangular matrix
Func11 <- function(A,b){
  n <- length(A[,1])
  for(k in 1:n) {
    pivot <- k                                    #Choose a pivot for Gaussian Elimination
    for(i in k:n){
      if(abs(A[pivot,k]) < abs(A[i,k]))
        pivot <- i
    }
    if(abs(A[pivot,k]) == 0){
      val <- matrix(0,nrow=n,ncol=1)              #Set a criterion so as to declare the matrix as singular or non-singular
      break
    }
    else{
      if(k!=n){                                   #We don't need to interchange the last row or do any computation regarding it
        val <- matrix(1,nrow=n,ncol=1)            #Set the criterion same as done in the previous loop
        Interchange <- A[pivot,]
        A[pivot,] <- A[k,]
        A[k,] <- Interchange                      #Interchange the k^th row with the pivot row for A
        b_Interchange <- b[pivot,1]
        b[pivot,1] <- b[k,1]
        b[k,1] <- b_Interchange                   #Interchange the k^th element with the pivot element for B
        for(i in (k+1):n){
          factor_ik <- A[i,k]/A[k,k]
          for(j in k:n){
            A[i,j] <- A[i,j] - factor_ik*A[k,j]   #Create an upper triangular matrix
          }
          b[i,1] <- b[i,1] - factor_ik*b[k,1]     #Make changes in b vector according to A matrix
        }
      }
    }
  }
  data1 <- cbind(A,b,val)                         #Return a new matrix with the elements of A, b and val
  return(data1)
}
#Compute the output vector
Func12 <- function(A,b){
  n <- length(A[,1])
  x <- matrix(1,nrow=n,ncol=1)
  x[n,1] <- b[n,1]/A[n,n]                         #Get the value of x_n
  for(i in (n-1):1){
    Sum <- 0
    for(j in (i+1):n){
      Sum <- Sum + A[i,j]*x[j,1]
    }
    x[i,1] <- (b[i,1] - Sum)/A[i,i]               #Get the values of x_i
  }
  return(x)
}
#Use the above defined function and taken the input as the coefficient matrix and the constant vector and give x as the output vector if the matrix is non-singular or give an error
Func13 <- function(Coeff,Constant){
  A <- Matrix_Conv_A(Coeff)
  b <- Matrix_Conv_B(Constant)
  c <- length(A[1,])
  data2 <- Func11(A,b)                         #Create an upper diagonal Matrix
  if(data2[1,c+2] == 0){
    print("Singular Matrix")
  } else{
    data2 <- data.frame(data2)
    New_A <- data2[,1:c]
    New_B <- data2[,(c+1)]
    New_B <- data.frame(New_B)
    print(Func12(New_A,New_B))
  }
}

########################################################
#Gauss-Jordan Method
#Define Functions
########################################################
#Create a upper triangular matrix
Func1 <- function(A,b){
  n <- length(A[,1])
  for(k in 1:n) {
    pivot <- k                                    #Choose a pivot for Gaussian Elimination
    for(i in k:n){
      if(abs(A[pivot,k]) < abs(A[i,k]))
        pivot <- i
    }
    if(abs(A[pivot,k]) == 0){
      val <- matrix(0,nrow=n,ncol=1)              #Set a criterion so as to declare the matrix as singular or non-singular
      break
    }
    else{
      if(k!=n){                                   #We don't need to interchange the last row or do any computation regarding it
        val <- matrix(1,nrow=n,ncol=1)            #Set the criterion same as done in the previous loop
        Interchange <- A[pivot,]
        A[pivot,] <- A[k,]
        A[k,] <- Interchange                      #Interchange the k^th row with the pivot row for A
        b_Interchange <- b[pivot,1]
        b[pivot,1] <- b[k,1]
        b[k,1] <- b_Interchange                   #Interchange the k^th element with the pivot element for B
        for(i in (k+1):n){
          factor_ik <- A[i,k]/A[k,k]
          for(j in k:n){
            A[i,j] <- A[i,j] - factor_ik*A[k,j]   #Create an upper triangular matrix
          }
          b[i,1] <- b[i,1] - factor_ik*b[k,1]     #Make changes in b vector according to A matrix
        }
      }
    }
  }
  data1 <- cbind(A,b,val)                         #Return a new matrix with the elements of A, b and val
  return(data1)
}
#Create a lower triangular matrix
Func4 <- function(A,b){
  n <- length(A[,1])
  for(k in n:2){
    for(i in (k-1):1){
      factor_ik <- A[i,k]/A[k,k]
      for(j in k:1){
        A[i,j] <- A[i,j] - factor_ik*A[k,j]             #Form the lower diagonal matrix
      }
      b[i,1] <- b[i,1] - factor_ik*b[k,1]               #Make corresponding changes in the outcome vector(b)
    }
  }
  data1 <- cbind(A,b)
  return(data1)
}
#Compute the output vector
Func5 <- function(A,b){
  n <- length(A[1,])
  x <- matrix(1,nrow=n,ncol=1)
  for(i in 1:n){
    x[i] <- (b[i,1])/A[i,i]                         #Compute the values of the vector x
  }
  return(x)
}
#Use the above defined function and taken the input as the coefficient matrix and the constant vector and give x as the output vector if the matrix is non-singular or give an error
Func3 <- function(Coeff,Constant){
  A <- Matrix_Conv_A(Coeff)
  b <- Matrix_Conv_B(Constant)
  c <- length(A[1,])
  data2 <- Func1(A,b)                         #Create an upper diagonal Matrix
  if(data2[1,c+2] == 0){
    print("Singular Matrix")
  } else{
    data2 <- data.frame(data2)
    New_A <- data2[,1:c]
    New_B <- data2[,(c+1)]
    New_B <- data.frame(New_B)
    data3 <- Func4(New_A,New_B)               #Create a lower diagonal matrix
    data3 <- data.frame(data3)
    A_New <- data3[,1:c]
    B_New <- data3[,(c+1)]
    B_New <- data.frame(B_New)
    print(Func5(A_New,B_New))                 #Compute the values of x if it exists
  }
}
######################################################
Func6 <- function(A,b,option){
  From_Package <- optR(Matrix_Conv_A(A),Matrix_Conv_B(b), method="LU")
  Value_Package <- From_Package$beta[1]
  Value_Package <- Value_Package[[1]]
  if(option == 1){
    Kase <- Func13(A,b)                                             
    print("Solved Using Gauss Elimination")               #Gauss-Elimination Method
  } else {
    Kase <- Func3(A,b)
    print("Solved using Gauss Jordan")                    #Gauss-Jordan Method  
  }
  error <- norm(Kase - Value_Package)/norm(Value_Package)
  cat("Relative error is ", error)
}
#For Code Coverage
Func7 <- function(A,b,option){
  if(option == 1){
    Func13(A,b)
  } else {
    Func3(A,b)
  }
}
######################################################
#Give the input vector provided by the user
A <- c(1,2,3,4,5,6,17,9,10,21,12,13,18,10,19,21)            #Give any coefficient matrix
b <- c(12,19,1,67)                                        #Give a constant vector
option <- 2
Func6(A,b,option)                                         #Get the value of x if it's singular or else get an error that it is singular
#####################################################