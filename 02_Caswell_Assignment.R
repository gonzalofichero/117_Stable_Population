#### 1.0. Mise en Place ####

##### Loading packages ##### 
library(tidyverse)
library(R.matlab)


#####  Importing data ##### 
pen1 <- as.matrix(read.table("penguinmat1.txt"))
pen2 <- as.matrix(read.table("penguinmat2.txt"))

# Checking
View(pen1)
View(pen2)


##### 1.a. Write a program to project the population starting from several different initial conditions: ##### 

# A) one newborn baby penguin
# B) one breeding adult penguin
# C) a population with the stable stage distribution


###### Extracting eigenvectors/values from A ###### 
eigen.penguin1 <- eigen(pen1)

###### Let's find the one that checks Perron-Frobenius ###### 
lambda1.position.pen1 <- which.max(abs(eigen.penguin1$values))

###### Eigenvalue + Rigth Eigenvector ###### 
lambda1.pen1 <- eigen.penguin1$values[lambda1.position.pen1]
omega1.pen1 <- eigen.penguin1$vectors[,lambda1.position.pen1]

###### LeftEigenvector ###### 
ve1 <- eigen(t(pen1))$vectors[,lambda1.position.pen1]


###### U matrix ###### 
pen1[1,7] <- 0
U.pen1 <- pen1

pen2[1,7] <- 0
U.pen2 <- pen2

###### Importing data ###### 
pen1 <- as.matrix(read.table("penguinmat1.txt"))
pen2 <- as.matrix(read.table("penguinmat2.txt"))


###### Calculating Fundamental Matrix N for both penguin Universes ###### 

N.pen1 <- solve(diag(dim(U.pen1)[1]) - U.pen1)

N.pen2 <- solve(diag(dim(U.pen2)[1]) - U.pen2)




#### 3. Population projection ####

##### Setting initial population in t = 0 #####
n0.a <- as.vector(c(1,0,0,0,0,0,0))
n0.b <- as.vector(c(0,0,0,0,0,0,1))
n0.c <- omega1.pen1


##### Calculating C vector #####
c1.a <- as.vector(ve1 %*% n0.a)
c1.b <- as.vector(ve1 %*% n0.b)
c1.c <- as.vector(ve1 %*% n0.c)



##### Looping for 10 periods #####
n.a <- matrix(,nrow=7,ncol=11)

for (t in 0:10) {
  n.a[,t+1] <- c1.a %*% (lambda1^t) %*% omega1
}

n.a <- as.data.frame(n.a)
names(n.complete) <- c("n0", "n1", "n2", "n3", "n4", "n5",
                       "n6", "n7", "n8", "n9", "n10")