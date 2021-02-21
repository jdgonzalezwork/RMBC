# RMBC package: Robust and efficient Clustering 
====================================================

#### Juan D. Gonzalez, Ricardo Maronna, Victor J. Yohai and Ruben H. Zamar 


## Introduction
------------

This package implements a clustering algorithm similar to Expectation Maximization 
algorithm for multivariate Gaussian finite mixture models by using Robust S estimator its main advantage is that The estimator is resistant to outliers, that means that results of
    estimator are still correct when there are atipycal values in the   sample.

The reference is _Robust Model-Based_ Clustering, Juan D. Gonzalez, Ricardo Maronna
, Victor J. Yohai, and Ruben H. Zamar [arxiv:2102.06851](https://arxiv.org/pdf/2102.06851.pdf).



# How to use the package RMBC
----------------------------------

``` {.r}
# Generate Sintetic data (three normal cluster in two dimension)
# clusters have different shapes and orentation.
# The data is contaminated uniformly (level 20%).
################################################
#### Start data generating process ############
##############################################

# generates base clusters

Z1 <- c(rnorm(100,0),rnorm(100,0),rnorm(100,0))
Z2 <- rnorm(300);
X <-  matrix(0, ncol=2,nrow=300);
X[,1]=Z1;X[,2]=Z2
true.cluster= c(rep(1,100),rep(2,100),rep(3,100))
# rotate, expand and translate base clusters
theta=pi/3;
aux1=matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow=2)

aux2=sqrt(4)*diag(c(1,1/4))

B=aux1%*%aux2%*%t(aux1)

X[true.cluster==3,]=X[true.cluster==3,]%*%aux2%*%aux1 + 
matrix(c(15,2), byrow = TRUE,nrow=100,ncol=2)
X[true.cluster==2,2] = X[true.cluster==2,2]*4
X[true.cluster==1,2] = X[true.cluster==1,2]*0.1
X[true.cluster==1, ] = X[true.cluster==1,]+ 
matrix(c(-15,-1),byrow = TRUE,nrow=100,ncol=2)

### Generate 60 sintetic outliers (contamination level 20%)

outliers=sample(1:300,60)
X[outliers, ] <- matrix(runif( 40, 2 * min(X), 2 * max(X) ),
                        ncol = 2, nrow = 60)

###############################################
#### END data generating process ############
#############################################

### APLYING RMBC ALGORITHM 

ret = RMBC(Y=X, K=3,max_iter = 82)

cluster = ret$cluster
#############################################
### plotting results ########################
#############################################
oldpar=par(mfrow=c(1,2))
plot(X,  main="actual clusters" )
for (j in 1:3){
  points(X[true.cluster==j,],pch=19, col=j+1)
}
points(X[outliers,],pch=19,col=1)

plot(X,main="clusters estimation")
for (j in 1:3){
  points(X[cluster==j,],pch=19, col=j+1)
}
points(X[ret$outliers,],pch=19,col=1)
par(oldpar)
```

The preprint [arxiv:1906.08198](https://arxiv.org/abs/1906.08198)) contains comparison with other robust model-based clustering procedures as well as technical details and applications.   


