# RMBC package: Robust Model-Based Clustering
====================================================

#### Juan D. Gonzalez, Ricardo Maronna, Victor J. Yohai and Ruben H. Zamar


## Introduction
------------

This package implements a clustering algorithm similar to Expectation Maximization
algorithm for multivariate Gaussian finite mixture models by using Robust S estimator.
Its main advantage is that The estimator is resistant to outliers, that means that results of
estimator are still correct when there are atipycal values in the   sample.

The reference is _Robust Model-Based_ Clustering, Juan D. Gonzalez, Ricardo Maronna
, Victor J. Yohai, and Ruben H. Zamar [arxiv:2102.06851](https://arxiv.org/pdf/2102.06851v2.pdf).



# How to use the package RMBC
------------------------------

Firstly, you should install the package from github, by typing

``` {.r}
devtools::install_github("jdgonzalezwork/RMBC")
```
where the package [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) must be previously installed.

The following routine generates synthetic data and runs the `RMBC` algorithm.

``` {.r}
library("RMBC")
# Generate Sintetic data (three normal cluster in two dimension)
# clusters have different shapes and orentation.
# The data is contaminated uniformly (level 20%).
################################################
#### Start data generating process ############
##############################################

# generates base clusters
set.seed(1)
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
X[true.cluster==2,2] = X[true.cluster==2,2]
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
plot(X,  main="Actual Clusters")
for (j in 1:3){
  points(X[true.cluster==j,],pch=19, col=j+1)
}
points(X[outliers,],pch=19,col=1)

plot(X,main="Clusters Estimation")
for (j in 1:3){
  points(X[cluster==j,],pch=19, col=j+1)
}
points(X[ret$outliers,],pch=19,col=1)
par(oldpar)
```
The graphical result of this routine can be viewed below, on the left image, the actual 3 clusters are in different colors, whereas the atypical values are represented in black. The right image shows how the algorithm was able to detect the three clusters and outliers.
Of course this is a small example, but in the reference [1], this routine is applied to
several synthetic scenarios as well as to real data.

![alt text](https://github.com/jdgonzalezwork/RMBC/blob/main/result.png)



The preprint [1] contains theoretical results and comparisons with other robust model-based clustering procedures. It also has technical details and applications to real data.

*References:*

[[1] Gonzalez, J. D., Maronna, R., Yohai, V. J., & Zamar, R. H. (2021). Robust Model-Based Clustering. arXiv preprint arXiv:2102.06851.](https://arxiv.org/pdf/2102.06851v2.pdf)
