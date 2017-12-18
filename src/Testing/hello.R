Funct <- function(A,b,option){
  n <- (length(A))^0.5
  A <- matrix(A,nrow = n,ncol=n,byrow = TRUE)
  n <- length(b)
  b <- matrix(b,nrow=n,ncol=1,byrow = TRUE)
  n <- length(A[,1])


  if(option==1){
    c <- length(A[1,])
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
        if(k!=n){                                   	#We don't need to interchange the last row or do any computation regarding it
          val <- matrix(1,nrow=n,ncol=1)            	#Set the criterion same as done in the previous loop
          Interchange <- A[pivot,]
          A[pivot,] <- A[k,]
          A[k,] <- Interchange                      	#Interchange the k^th row with the pivot row for A
          b_Interchange <- b[pivot,1]
          b[pivot,1] <- b[k,1]
          b[k,1] <- b_Interchange                   	#Interchange the k^th element with the pivot element for B
          for(i in (k+1):n){
            factor_ik <- A[i,k]/A[k,k]
            for(j in k:n){
              A[i,j] <- A[i,j] - factor_ik*A[k,j]   	#Create an upper triangular matrix
            }
            b[i,1] <- b[i,1] - factor_ik*b[k,1]     	#Make changes in b vector according to A matrix
          }
        }
      }
    }
    data2 <- cbind(A,b,val)                         	#Return a new matrix with the elements of A, b and val
    if(data2[1,c+2] == 0){
      print("Singular Matrix")
    } else {
      data2 <- data.frame(data2)
      A <- data2[,1:c]
      b <- data2[,(c+1)]
      b <- data.frame(b)
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
  }
  if(option == 2){
    c <- length(A[1,])
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
    data2 <- cbind(A,b,val)                         #Return a new matrix with the elements of A, b and val
    if(data2[1,c+2] == 0){
      print("Singular Matrix")
    } else {
      data2 <- data.frame(data2)
      A <- data2[,1:c]
      b <- data2[,(c+1)]
      b <- data.frame(b)
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
      data3 <- cbind(A,b)
      A <- data3[,1:c]
      b <- data3[,(c+1)]
      b <- data.frame(b)
      n <- length(A[1,])
      x <- matrix(1,nrow=n,ncol=1)
      for(i in 1:n){
        x[i] <- (b[i,1])/A[i,i]                         #Compute the values of the vector x
      }
      return(x)
    }
  }
}
