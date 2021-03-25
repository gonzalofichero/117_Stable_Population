#### 0.0. Mise en Place ####

##### Loading packages ##### 
library(tidyverse)
library(R.matlab)


#####  Importing data ##### 
pen1 <- as.matrix(read.table("penguinmat1.txt"))
pen2 <- as.matrix(read.table("penguinmat2.txt"))

# Checking
View(pen1)
View(pen2)


#### 1.a. Write a program to project the population starting from several different initial conditions: #### 

# A) one newborn baby penguin
# B) one breeding adult penguin
# C) a population with the stable stage distribution


##### Normal Population: penguinmat1.txt #####

###### Extracting eigenvectors/values from A ###### 
eigen.penguin1 <- eigen(pen1)

###### Let's find the one that checks Perron-Frobenius ###### 
lambda1.position.pen1 <- which.max(abs(eigen.penguin1$values))

###### Eigenvalue + Rigth Eigenvector ###### 
lambda1.pen1 <- eigen.penguin1$values[lambda1.position.pen1]
omega1.pen1 <- eigen.penguin1$vectors[,lambda1.position.pen1]

# Standardizing omega vector
omega1.pen1.stand <- as.numeric(omega1.pen1) / as.numeric(sum(omega1.pen1))

###### LeftEigenvector ###### 
ve1 <- eigen(t(pen1))$vectors[,lambda1.position.pen1]


###### U matrix ###### 
U.pen1 <- pen1
U.pen1[1,7] <- 0


###### Calculating Fundamental Matrix N for both penguin Universes ###### 

N.pen1 <- solve(diag(dim(U.pen1)[1]) - U.pen1)


###### Population projection ###### 

####### Setting initial population in t = 0 ####### 
n0.a <- as.vector(c(1,0,0,0,0,0,0))
n0.b <- as.vector(c(0,0,0,0,0,0,1))
n0.c <- omega1.pen1.stand


####### Calculating C vector ####### 
c1.a <- as.vector(ve1 %*% n0.a)
c1.b <- as.vector(ve1 %*% n0.b)
c1.c <- as.vector(ve1 %*% n0.c)


####### Looping for X periods ####### 
periods <- 25
n.a <- matrix(,nrow=dim(U.pen1)[1],ncol=periods+1)
n.b <- matrix(,nrow=dim(U.pen1)[1],ncol=periods+1)
n.c <- matrix(,nrow=dim(U.pen1)[1],ncol=periods+1)

####### Adding initial population #######
n.a[,1] <- n0.a
n.b[,1] <- n0.b
n.c[,1] <- n0.c

####### Looping #######
for (t in 1:periods) {
  
  n.a[,t+1] <- pen1 %*% n.a[,t]
  n.b[,t+1] <- pen1 %*% n.b[,t]
  n.c[,t+1] <- pen1 %*% n.c[,t]
  
}


###### Re-arranging projections ######

# Initial population a #
n.a.dat <- data.frame(n.a)
names(n.a.dat) <- seq(0,periods)
n.a.dat$stage <- seq(1,dim(U.pen1)[1])

n.a.dat <- n.a.dat %>% 
            pivot_longer(cols = !stage, values_to = "n_t", names_to = "t") %>% 
            mutate(stage = as.factor(stage),
                   t = as.numeric(t),
                   type_penguin = "normal",
                   initial_pop = "a")

# Initial population b #
n.b.dat <- data.frame(n.b)
names(n.b.dat) <- seq(0,periods)
n.b.dat$stage <- seq(1,dim(U.pen1)[1])

n.b.dat <- n.b.dat %>% 
  pivot_longer(cols = !stage, values_to = "n_t", names_to = "t") %>% 
  mutate(stage = as.factor(stage),
         t = as.numeric(t),
         type_penguin = "normal",
         initial_pop = "b")

# Initial population c #
n.c.dat <- data.frame(n.c)
names(n.c.dat) <- seq(0,periods)
n.c.dat$stage <- seq(1,dim(U.pen1)[1])

n.c.dat <- n.c.dat %>% 
  pivot_longer(cols = !stage, values_to = "n_t", names_to = "t") %>% 
  mutate(stage = as.factor(stage),
         t = as.numeric(t),
         type_penguin = "normal",
         initial_pop = "c")


# Plotting all together #

penguins <- rbind(n.a.dat, n.b.dat, n.c.dat)

penguins %>% 
  ggplot(aes(x = t, y = n_t, color = stage)) + geom_line() +
  facet_wrap(~ initial_pop)










