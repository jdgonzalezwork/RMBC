  #'
  #' Robust Model Base Clustering  a robust and efficient version of EM algorithm.
  #' @param Y A matrix of size n x p.
  #' @param K The number of clusters.
  #' @param max_iter a maximum number of iterations used for the
  #' algorithm stopping rule
  #' @param tolerance tolerance parameter used for the algorithm stopping
  #'     rule
  #' @return A list including the estimated mixture distribution parameters and cluster-label for the
  #' observations  
  #' \itemize{
  #' \item{\code{alpha}}{: K numeric values representing the convex combination coefficients.}
  #' \item{\code{mu}}{: a list of length K with the location initial estimators.}
  #' \item{\code{sigma}}{: a list of length K with the location scatter matrix estimators.}
  #' \item{\code{nonoutliers}}{: an array of indices that contains the estimated nonoutliers observations}
  #' \item{\code{outliers}}{: an array of indices that contains the estimated outliers observations}
  #' }
   
  #' @examples 
  #' # Generate Sintetic data (three normal cluster in two dimension)
  #' # clusters have different shapes and orentation.
  #' # The data is contaminated uniformly (level 20%).
  #' ################################################
  #' #### Start data generating process ############
  #' ##############################################
  #' 
  #' # generates base clusters
  #' 
  #' Z1 <- c(rnorm(100,0),rnorm(100,0),rnorm(100,0))
  #' Z2 <- rnorm(300);
  #' X <-  matrix(0, ncol=2,nrow=300);
  #' X[,1]=Z1;X[,2]=Z2
  #' true.cluster= c(rep(1,100),rep(2,100),rep(3,100))
  #' # rotate, expand and translate base clusters
  #' theta=pi/3;
  #' aux1=matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow=2)
  #' 
  #' aux2=sqrt(4)*diag(c(1,1/4))
  #' 
  #' B=aux1%*%aux2%*%t(aux1)
  #' 
  #' X[true.cluster==3,]=X[true.cluster==3,]%*%aux2%*%aux1 + 
  #' matrix(c(15,2), byrow = TRUE,nrow=100,ncol=2)
  #' X[true.cluster==2,2] = X[true.cluster==2,2]*4
  #' X[true.cluster==1,2] = X[true.cluster==1,2]*0.1
  #' X[true.cluster==1, ] = X[true.cluster==1,]+ 
  #' matrix(c(-15,-1),byrow = TRUE,nrow=100,ncol=2)
  #' 
  #' ### Generate 60 sintetic outliers (contamination level 20%)
  #' 
  #' outliers=sample(1:300,60)
  #' X[outliers, ] <- matrix(runif( 40, 2 * min(X), 2 * max(X) ),
  #'                         ncol = 2, nrow = 60)
  #' 
  #' ###############################################
  #' #### END data generating process ############
  #' #############################################
  #' 
  #' ### APLYING RMBC ALGORITHM 
  #' 
  #' ret = RMBC(Y=X, K=3,max_iter = 82)
  #' 
  #' cluster = ret$cluster
  #' #############################################
  #' ### plotting results ########################
  #' #############################################
  #' oldpar=par(mfrow=c(1,2))
  #' plot(X,  main="actual clusters" )
  #' for (j in 1:3){
  #'   points(X[true.cluster==j,],pch=19, col=j+1)
  #' }
  #' points(X[outliers,],pch=19,col=1)
  #' 
  #' plot(X,main="clusters estimation")
  #' for (j in 1:3){
  #'   points(X[cluster==j,],pch=19, col=j+1)
  #' }
  #' points(X[ret$outliers,],pch=19,col=1)
  #' par(oldpar)
  #'
  #' @export
  
  RMBC=function(Y,K,max_iter=80, tolerance=1e-4){
  
    niterFixedPoint=1;
    result=robustINIT(Y = Y,K = K)
    resultRMBCaux=RMBCaux(Y = Y,K = K,
                          thetaOld.alpha = result$alphaINIT,
                          thetaOld.mu=result$muINIT,
                          thetaOld.sigma = result$sigmaINIT,
                          max_iter,niterFixedPoint,tolerance)
    
    Sigma=resultRMBCaux$theta.sigma;
    mu=resultRMBCaux$theta.mu;
    alpha=resultRMBCaux$theta.alpha
    outliers=resultRMBCaux$outliers
    nonoutliers=resultRMBCaux$nonoutliers
    iter=resultRMBCaux$iter
    cluster=resultRMBCaux$cluster
    
    list(mu=mu,Sigma=Sigma,alpha=alpha, outliers=outliers,
         nonoutliers=nonoutliers, cluster=cluster,iter=iter)
  }
  