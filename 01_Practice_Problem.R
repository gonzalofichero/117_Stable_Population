#### 0. Mise en Place ####

##### Loading packages ##### 
library(tidyverse)
library(R.matlab)


#####  Importing Matlab data into R ##### 
population <- readMat("dicymbe_mat.mat")

#####  Checking: ##### 
glimpse(population)

# Data is an object with F, U and A matrices
# Let's take each apart

##### Extracting matrices ##### 
A <- as.matrix(population$a)
F <- as.matrix(population$f)
U <- as.matrix(population$u)


#### 1. Extracting eigenvectors/values from A ####

eigen.A <- eigen(A)

eigen.A$values

#####   Let's find the one that checks Perron-Frobenius #####  
lambda1.position <- which.max(abs(eigen.A$values))

#####   Eigenvalue + Rigth Eigenvector #####  
lambda1 <- eigen.A$values[lambda1.position]
omega1 <- (-1) * eigen.A$vectors[,lambda1.position]

##### LeftEigenvector #####  
ve1 <- eigen(t(A))$vectors[,lambda1.position]



#### 2. Calculating Fundamental Matrix N ####

N <- solve(diag(15) - U)



#### 3. Population projection ####

##### Setting initial population in t = 0 #####
n0 <- as.vector(rep(10000, 15))

##### Calculating C vector #####
c1 <- as.vector(ve1 %*% n0)


##### Looping for 10 periods #####
n <- matrix(,nrow=15,ncol=11)

for (t in 0:10) {
  n[,t+1] <- c1 %*% (lambda1^t) %*% omega1
}

n.complete <- as.data.frame(n)
names(n.complete) <- c("n0", "n1", "n2", "n3", "n4", "n5",
                       "n6", "n7", "n8", "n9", "n10")




