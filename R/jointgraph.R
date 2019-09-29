#' @title Joint Estimation of Multiple Graphical Models

#' @description The method in paper Joint Estimation of Multiple Graphical Models
#' A method which jointly estimates several graphical models corresponding to the different categories
#' present in the data. The method aims to preserve the common structure, while allowing for
#' differences between the categories.
#'
#' @param trainX: data
#' @param trainY: labels for categories (1, 2, 3,...)
#' @param lambda_value: tuning parameter
#' @param adaptive:
#'
#' @return outputs: parameters of the model trained with data trainX and trainY

#' @examples
#' print("this is a test for CGM_AHP_train")
#' X = matrix(rnorm(50*20),ncol=20)
#' Y <- NULL
#' ps <- cumsum(c(0.25,0.25,0.25,0.25))
#' r <- runif(50)
#' x <- c(1,2,3,4)
#' for (i in 1:50) Y <- c(Y, x[which(r[i] <= ps)[1]])
#' trainX = X
#' trainY = Y
#' lbd = 0.001
#' result = jointgraph_train(trainX, trainY, lbd)
#'
#' @export jointgraph_train

library(glasso)

jointgraph_train <- function(trainX, trainY, lambda_value, adaptive_weight=array(1, c(length(unique(Y)), ncol(X), ncol(X))))
{
  ## Set the general paramters
  K <- length(unique(trainY))
  p <- ncol(trainX)
  diff_value <- 1e+10
  count <- 0
  tol_value <- 1e-2
  max_iter <- 30

  ## Set the optimizaiton parameters
  OMEGA <- array(0, c(K, p, p))
  S <- array(0, c(K, p, p))
  OMEGA_new <- array(0, c(K, p, p))
  nk <- rep(0, K)

  ## Initialize Omega
  for (k in seq(1, K))
  {
    idx <- which(trainY == k)
    S[k, , ] <- cov(trainX[idx, ])
    if (kappa(S[k, , ]) > 1e+15)
    {
      S[k, , ] <- S[k, , ] + 0.001*diag(p)
    }
    tmp <- solve(S[k, , ])
    OMEGA[k, , ] <- tmp
    nk[k] <- length(idx)
  }

  ## Start loop
  while((count < max_iter) & (diff_value > tol_value))
  {
    tmp <- apply(abs(OMEGA), c(2,3), sum)
    tmp[abs(tmp) < 1e-10] <- 1e-10
    V <- 1 / sqrt(tmp)

    for (k in seq(1, K))
    {
      penalty_matrix <- lambda_value * adaptive_weight[k, , ] * V
      obj_glasso <- glasso::glasso(S[k, , ], penalty_matrix, maxit=100)
      OMEGA_new[k, , ] <- (obj_glasso$wi + t(obj_glasso$wi)) / 2
      #OMEGA_new[k, , ] <- obj_glasso$wi
    }

    ## Check the convergence
    diff_value <- sum(abs(OMEGA_new - OMEGA)) / sum(abs(OMEGA))
    count <- count + 1
    OMEGA <- OMEGA_new
    #cat(count, ', diff_value=', diff_value, '\n')
  }

  ## Filter the noise
  for (k in seq(1, K))
  {
    ome <- OMEGA[k, , ]
    ww <- diag(ome)
    ww[abs(ww) < 1e-10] <- 1e-10
    ww <- diag(1/sqrt(ww))
    tmp <- ww %*% ome %*% ww
    ome[abs(tmp) < 1e-5] <- 0
    OMEGA[k, , ] <- ome
  }

  output <- list()
  output$OMEGA <- OMEGA
  output$S <- S
  output$lambda <- lambda_value

  return(output)
}
