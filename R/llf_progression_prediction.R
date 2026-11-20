#' Predictions from a local linear progression forest
#'
#' Gets estimates of the conditional median of \eqn{Y} given \eqn{X} based on the trained progression forest.
#'
#' @param object The trained progression forest.
#' @param newdata A double vector, matrix or data.frame with test point where predictions should be made. It should have the the same number of columns as the training matrix \eqn{X}. Moreover, columns must appear in the same order.
#' @param pred.median A boolean indicating whether the conditional median should be computed or not. The default value is `TRUE`. Set to `FALSE` it uses the conditional mean loss function on standarized data to compute predictions before transforming them back to the original scale. This implementation will be faster in practice when \eqn{p > 1}, where \eqn{p} is the dimension of the dimension of a predictor \eqn{X}.
#' @param quantile A probability or a vector of probabilities with values in \eqn{(0,1)} where conditional quantiles should be computed.
#' @param lambda A double value with positive values computing a L1 penalization for quantile regression. Default value set to zero.
#' @param a A double with values in \eqn{[-1 , 1]} given the known value of the regression approximation \eqn{ $Y* \\approx a X*$  }
#' @param ... Additional arguments.
#'
#' @return A vector of prediction for the conditional median, or a matrix of quantile predictions.
#' @export
#'
#' @examples
#' # example code
#'
#' # set.seed(28)
#' n <- 1000
#' p <- 1
#' X <- matrix(rnorm(n * p), n, p)
#' reg_fun <- function(x) x^3/3
#' Y <- reg_fun(X[, 1]) + rnorm(n)
#' X_test <- matrix(rnorm(50 * p), 50, p)
#' X_test[,1] <- seq(-4,4, length.out=50)
#' # progression predictions
#' llf_progression_forest <- llf_progression(X, Y)
#' llf_predictions <- predict(llf_progression_forest, X_test)
#' curve(reg_fun, from = -4, to = 4)
#' points(X_test[,1],llf_predictions$predictions, lty=2)
#' rug(X[,1], side = 1, col="black")

predict.llf_progression <-
  function(object,
           newdata,
           pred.median = TRUE,
           quantile = NULL,
           lambda = 0,
           a = NULL,
           ...) {
    newdata <- as.matrix(newdata)
    p <- ncol(newdata)
    newdatatrans <- NULL

    both.sided <- object$both.sided
    old_upper <- object$upper
    upper <- NULL
    if (both.sided) {
      up <- (newdata[, 1] > object$upper_X)
      low <- (newdata[, 1] < object$lower_X)
      mid   <- !(up + low)
      upper <- list(upper = up,
                    lower = low,
                    mid = mid)
    }

    ## transforms new data to the Laplace scale
    for (i in 1:p)
      newdatatrans <-
      cbind(newdatatrans,
            object$trans_X[[i]]$data_to_laplace(newdata[, i], upper = upper))
    newdatatrans <- apply(newdatatrans, 2, bound_exp)
    Y_train_LP <- as.vector(object$Y_train_LP)
    X_train_LP <- as.matrix(object$X_train_LP)
    ## computes the weights from the local linear fit

    weix <- upper
    if (both.sided) {
      if(length((newdatatrans[upper$upper, ])) >0 ){
        weix$upper  <-
          grf::get_forest_weights(object$llf_LP$upper, newdata = as.matrix(newdatatrans[upper$upper, ]))
      }
      if(length((newdatatrans[upper$lower, ])) >0 ){
        weix$lower  <-
          grf::get_forest_weights(object$llf_LP$lower, newdata = as.matrix(newdatatrans[upper$lower, ]))
      }
      if(length((newdatatrans[upper$mid, ])) >0 ){
        weix$mid    <-
          grf::get_forest_weights(object$llf_LP$mid,  newdata = as.matrix(newdatatrans[upper$mid, ]))
      }
    } else{
      weix <-
        grf::get_forest_weights(object$llf_LP, newdata = newdatatrans)
    }

    if (!is.null(quantile) & all(quantile == .5)) {
      quantile = NULL
      pred.median = TRUE
    }
    if (pred.median == TRUE) {
      if (p > 1) {
        fn      <-
          function(par, x, wx)
            sum(wx * abs(Y_train_LP -  times(
              sapply(1:p, function(k)
                X_train_LP[, k] - x[k]) , 1 - 2 / (1 + exp(par[2:(p + 1)]))
            )  - par[1])) + lambda * sum(abs(1 - 2 / (1 + exp(par[2:(p + 1)]))))
        result  <-
          object$trans_Y$laplace_to_data(sapply(1:nrow(newdatatrans), function(k)
            stats::optim(rep(0, (p + 1)), function(par)
              fn(par, newdatatrans[k, ], weix[k, ]), method = "L-BFGS-B" )$par[1]), upper = upper)
      }
      else{
        if (both.sided) {
          if(!is.null(a)){
            graphics::par(mfrow=c(1,3))
            compute_results <- function(wei,
                                        newdt,
                                        Y_train_LP,
                                        X_train_LP,
                                        upp, a) {
              n        <- nrow(newdt)
              result   <- 1:n
              fn       <-
                function(par, x, wx)
                  sum(wx * abs(Y_train_LP -  a * (X_train_LP - x)  - par[1]) ^ 2 )
              init     <-
                sapply(1:n, function(k)
                  stats::optim(rep(0, (p + 1)), function(par)
                    fn(par, newdt[k, ], wei[k, ]), method = "L-BFGS-B" )$par)
              ##
              ind1    <- which(newdt[, 1] > 1 | newdt[, 1] < -1)
              ind2    <- which(!(1:n) %in% ind1)
              fn1     <-
                function(par, x, wx)
                  sum(wx * abs(Y_train_LP -  a * (X_train_LP - x)  - par[1]))
              fn2     <-
                function(par, x, wx)
                  sum(wx * abs(Y_train_LP -  a * (X_train_LP - x)  - par[1]))
              est1    <-
                sapply(ind1, function(k)
                  stats::optim(init[, k], function(par)
                    fn1(par, newdt[k, ], wei[k, ]), method = "L-BFGS-B")$par)
              est2    <-
                sapply(ind2, function(k)
                  stats::optim(init[, k], function(par)
                    fn2(par, newdt[k, ], wei[k, ]), method = "L-BFGS-B")$par)

              if (length(ind1) > 0)
                result[ind1]  <-
                object$trans_Y$laplace_to_data(est1[1, ],
                                               upper = list(
                                                 upper = upp$upper[ind1],
                                                 lower = upp$lower[ind1],
                                                 mid = upp$mid[ind1]
                                               ))
              if (length(ind2) > 0)
                result[ind2]  <-
                object$trans_Y$laplace_to_data(est2[1, ],
                                               upper = list(
                                                 upper = upp$upper[ind2],
                                                 lower = upp$lower[ind2],
                                                 mid = upp$mid[ind2]
                                               ))
              return(unlist(result))
            }
            n <- nrow(newdata)
            result     <- rep(0, n)

            ## upper
            Y_train_LP <- as.vector(object$Y_train_LP$upper)
            X_train_LP <- as.matrix(object$X_train_LP$upper)
            wei   <- as.matrix(weix$upper)
            newdt <- as.matrix(newdatatrans[upper$upper,])
            if (nrow(newdt) > 0) {
              upp   <-
                list(
                  upper = rep(T, length(newdt)),
                  lower = rep(F, length(newdt)),
                  mid = rep(F, length(newdt))
                )
              new   <-
                compute_results(wei, newdt, Y_train_LP, X_train_LP, upp, a = a)
              result[upper$upper]   <- new
            }

            ## plot for verification
            #curve(reg_fun, from = 0, to =  30)
            #points( newdata[upper$upper, ], result[upper$upper] )


            ## lower
            #plot(X_train_LP,Y_train_LP)
            #abline(0,1)
            Y_train_LP <- as.vector(object$Y_train_LP$lower)
            X_train_LP <- as.matrix(object$X_train_LP$lower)
            wei   <- as.matrix(weix$lower)
            newdt <- as.matrix(newdatatrans[upper$lower,])
            if (nrow(newdt) > 0) {
              upp   <-
                list(
                  upper = rep(F, length(newdt)),
                  lower = rep(T, length(newdt)),
                  mid = rep(F, length(newdt))
                )
              result[upper$lower]   <-
                compute_results(wei, newdt, Y_train_LP, X_train_LP, upp, a=a)
            }
            ## plot for verification
            #curve(reg_fun, from = -30, to =  0)
            #points( newdata[upper$lower, ], result[upper$lower] )

            ## mid
            Y_train_LP <- as.vector(object$Y_train_LP$mid)
            X_train_LP <- as.matrix(object$X_train_LP$mid)
            wei   <- as.matrix(weix$mid)
            newdt <- as.matrix(newdatatrans[upper$mid,])
            if (nrow(newdt) > 0) {
              upp   <-
                list(
                  upper = rep(F, length(newdt)),
                  lower = rep(F, length(newdt)),
                  mid = rep(T, length(newdt))
                )
              result[upper$mid]   <-
                compute_results(wei, newdt, Y_train_LP, X_train_LP, upp, a = 0)
            }
            ## plot for verification
            #curve(reg_fun, from = -5, to =  5)
            #points( newdata[upper$mid, ], result[upper$mid] )
            #plot(result[upper$mid],Y_train[upper$mid] )
            #abline(0,1)

          }else{
            compute_results_a <- function(wei,
                                          newdt,
                                          Y_train_LP,
                                          X_train_LP,
                                          upp) {
              n        <- nrow(newdt)
              result   <- 1:n
              fn       <-
                function(par, x, wx)
                  sum(wx * abs(Y_train_LP -  (1 - 2 * 1 / (
                    1 + exp(par[2])
                  )) * (X_train_LP - x)  - par[1])^ 2)
              init     <-
                sapply(1:n, function(k)
                  stats::optim(rep(0, (p + 1)), function(par)
                    fn(par, newdt[k, ], wei[k, ]), method = "L-BFGS-B")$par)
              ##
              ind1    <- which(newdt[, 1] > 1 | newdt[, 1] < -1)
              ind2    <- which(!(1:n) %in% ind1)
              fn1     <-
                function(par, x, wx)
                  sum(wx * abs(Y_train_LP -  (1 - 2 * 1 / (
                    1 + exp(par[2])
                  )) * (X_train_LP - x)  - par[1]))
              fn2     <-
                function(par, x, wx)
                  sum(wx * abs(Y_train_LP -  (1 - 2 * 1 / (
                    1 + exp(par[2])
                  )) * (X_train_LP - x)  - par[1]))
              est1    <-
                sapply(ind1, function(k)
                  stats::optim(init[, k], function(par)
                    fn1(par, newdt[k, ], wei[k, ]), method = "L-BFGS-B")$par)
              est2    <-
                sapply(ind2, function(k)
                  stats::optim(init[, k], function(par)
                    fn2(par, newdt[k, ], wei[k, ]), method = "L-BFGS-B")$par)

              if (length(ind1) > 0)
                result[ind1]  <-
                object$trans_Y$laplace_to_data(est1[1, ],
                                               upper = list(
                                                 upper = upp$upper[ind1],
                                                 lower = upp$lower[ind1],
                                                 mid = upp$mid[ind1]
                                               ))
              if (length(ind2) > 0)
                result[ind2]  <-
                object$trans_Y$laplace_to_data(est2[1, ],
                                               upper = list(
                                                 upper = upp$upper[ind2],
                                                 lower = upp$lower[ind2],
                                                 mid = upp$mid[ind2]
                                               ))
              return(unlist(result))
            }
            n <- nrow(newdata)
            result     <- rep(0, n)

            ## upper
            Y_train_LP <- as.vector(object$Y_train_LP$upper)
            X_train_LP <- as.matrix(object$X_train_LP$upper)
            wei   <- as.matrix(weix$upper)
            newdt <- as.matrix(newdatatrans[upper$upper,])
            if (nrow(newdt) > 0) {
              upp   <-
                list(
                  upper = rep(T, length(newdt)),
                  lower = rep(F, length(newdt)),
                  mid = rep(F, length(newdt))
                )
              new   <-
                compute_results_a(wei, newdt, Y_train_LP, X_train_LP, upp)
              result[upper$upper]   <- new
            }

            ## plot for verification
            #par(mfrow=c(1,2))
            #curve(reg_fun, from = 5, to =  30)
            #points( newdata[upper$upper, ], result[upper$upper] )
            #plot(result[upper$upper],Y_train[upper$upper] )
            #abline(0,1)
            ## lower
            #plot(X_train_LP,Y_train_LP)
            Y_train_LP <- as.vector(object$Y_train_LP$lower)
            X_train_LP <- as.matrix(object$X_train_LP$lower)
            wei   <- as.matrix(weix$lower)
            newdt <- as.matrix(newdatatrans[upper$lower,])
            if (nrow(newdt) > 0) {
              upp   <-
                list(
                  upper = rep(F, length(newdt)),
                  lower = rep(T, length(newdt)),
                  mid = rep(F, length(newdt))
                )
              result[upper$lower]   <-
                compute_results_a(wei, newdt, Y_train_LP, X_train_LP, upp)
            }
            ## plot for verification
            #curve(reg_fun, from = -30, to =  -5)
            #points( newdata[upper$lower, ], result[upper$lower] )
            #plot(result[upper$lower],Y_train[upper$lower] )
            #abline(0,1)
            ## mid
            Y_train_LP <- as.vector(object$Y_train_LP$mid)
            X_train_LP <- as.matrix(object$X_train_LP$mid)
            wei   <- as.matrix(weix$mid)
            newdt <- as.matrix(newdatatrans[upper$mid,])
            if (nrow(newdt) > 0) {
              upp   <-
                list(
                  upper = rep(F, length(newdt)),
                  lower = rep(F, length(newdt)),
                  mid = rep(T, length(newdt))
                )
              result[upper$mid]   <-
                compute_results_a(wei, newdt, Y_train_LP, X_train_LP, upp)
            }
            ## plot for verification
            #curve(reg_fun, from = -15, to =  15)
            #points( newdata[upper$mid, ], result[upper$mid] )
            #plot(result[upper$mid],Y_train[upper$mid] )
            #abline(0,1)

          }

        } else{
          n        <- nrow(newdata)
          fn       <-
            function(par, x, wx)
              sum(wx * abs(Y_train_LP -  (1 - 2 * 1 / (
                1 + exp(par[2])
              )) * (X_train_LP - x)  - par[1]) ^ 2)
          init     <-
            sapply(1:n, function(k)
              stats::optim(rep(0, (p + 1)), function(par)
                fn(par, newdatatrans[k, ], weix[k, ]), method = "L-BFGS-B")$par)

          ## fast
          ind1    <-
            which(newdatatrans[, 1] > 1.2 | newdatatrans[, 1] < -1.2)
          ind2    <- which(!(1:n) %in% ind1)
          fn1     <-
            function(par, x, wx)
              sum(wx * abs(Y_train_LP -  (1 - 2 * 1 / (
                1 + exp(par[2])
              )) * (X_train_LP - x)  - par[1]))
          fn2     <-
            function(par, x, wx)
              sum(wx * abs(Y_train_LP -  (1 - 2 * 1 / (
                1 + exp(par[2])
              )) * (X_train_LP - x)  - par[1]))
          est1    <-
            sapply(ind1, function(k)
              stats::optim(init[, k], function(par)
                fn1(par, newdatatrans[k, ], weix[k, ]), method = "L-BFGS-B")$par)
          est2    <-
            sapply(ind2, function(k)
              stats::optim(init[, k], function(par)
                fn2(par, newdatatrans[k, ], weix[k, ]), method = "L-BFGS-B")$par)
          result     <- 1:n

          if (length(ind1) > 0)
            result[ind1]  <-
            object$trans_Y$laplace_to_data(est1[1, ], upper = upper)
          if (length(ind2) > 0)
            result[ind2]  <-
            object$trans_Y$laplace_to_data(est2[1, ], upper = upper)
        }

      }
    } else if (pred.median == FALSE) {
      if (both.sided) {
        compute_results_quantile <- function(wei,
                                             newdt,
                                             Y_train_LP,
                                             X_train_LP,
                                             upp) {
          n        <- nrow(newdt)
          result   <- 1:n

          if (is.null(quantile))
            fn       <-
            function(par, x, wx)
              sum(wx * abs(Y_train_LP -  (1 - 2 * 1 / (
                1 + exp(par[2])
              )) * (X_train_LP - x)  - par[1]) ^ 2)
          else
            fn       <-
            function(par, x, wx)
              sum(wx * quantile_loss(Y_train_LP ,  (1 - 2 * 1 / (
                1 + exp(par[2])
              )) * (X_train_LP - x)  - par[1], tau=quantile))
          init     <-
            sapply(1:n, function(k)
              stats::optim(rep(0, (p + 1)), function(par)
                fn(par, newdt[k, ], wei[k, ]), method = "L-BFGS-B")$par)
          result   <- object$trans_Y$laplace_to_data(init[1, ], upp)
          return(unlist(result))
        }
        n <- length(newdata)
        result     <- rep(0, n)

        ## upper
        Y_train_LP <- as.vector(object$Y_train_LP$upper)
        X_train_LP <- as.matrix(object$X_train_LP$upper)
        wei   <- as.matrix(weix$upper)
        newdt <- as.matrix(newdatatrans[upper$upper,])
        upp   <-
          list(
            upper = rep(T, length(newdt)),
            lower = rep(F, length(newdt)),
            mid = rep(F, length(newdt))
          )
        new   <-
          compute_results_quantile(wei, newdt, Y_train_LP, X_train_LP, upp)
        result[upper$upper]   <- new

        ## lower
        #plot(X_train_LP,Y_train_LP)
        Y_train_LP <- as.vector(object$Y_train_LP$lower)
        X_train_LP <- as.matrix(object$X_train_LP$lower)
        wei   <- as.matrix(weix$lower)
        newdt <- as.matrix(newdatatrans[upper$lower,])
        upp   <-
          list(
            upper = rep(F, length(newdt)),
            lower = rep(T, length(newdt)),
            mid = rep(F, length(newdt))
          )
        result[upper$lower]   <-
          compute_results_quantile(wei, newdt, Y_train_LP, X_train_LP, upp)

        ## mid
        Y_train_LP <- as.vector(object$Y_train_LP$mid)
        X_train_LP <- as.matrix(object$X_train_LP$mid)
        wei   <- as.matrix(weix$mid)
        newdt <- as.matrix(newdatatrans[upper$mid,])
        upp   <-
          list(
            upper = rep(F, length(newdt)),
            lower = rep(F, length(newdt)),
            mid = rep(T, length(newdt))
          )
        result[upper$mid]   <-
          compute_results_quantile(wei, newdt, Y_train_LP, X_train_LP, upp)

      } else{
        pred    <-
          utils::getFromNamespace("predict.ll_regression_forest", "grf")
        result  <-
          object$trans_Y$laplace_to_data(bound_exp(pred(object$llf_LP, newdatatrans)$predictions), upper =
                                           upper)
      }

    }
    return(list(predictions = result, quantile = quantile))
  }



#' Bounds for the Laplace transformation
#'
#' Bound the entries of each vector by \eqn{-36} and \eqn{36} as for \eqn{x} above or below these bound \eqn{exp(36) = Inf}.
#'
#' @param x double vector to be bounded.
#'
#' @return vector with bounded entries
#'
#' @examples
#' #x <- seq(-100,100,by=5)
#' #xp <- bound_exp(x)
bound_exp  <- function(x){
  sapply(x, function(k) max(-36, min(k,36) ))
}
#' Column-wise parameter multiplication and addition
#'
#' Computes a dot product between each matrix row and the vector of parameters
#'
#' @param mat A double matrix of dimension n x p
#' @param par A vector of parameters of dimension p
#'
#' @return An vector of dimension c
#'
#' @examples
#' #n <- 10
#' #p <- 3
#' #par <- c(1,-3,.5)
#' #X <- matrix(rnorm(n * p), ncol = p )
#' #res <- times( X, par )
times      <- function(mat,par){
  rowSums( mat %*% par )
}

#' Computes the upper and lower bound restrictions for optimization
#'
#'@description
#'  `compute_ui()` returns the matrix of upper bounds
#'  `compute_ci()` returns the matrix of lower bounds
#' @param p An integer vector indicating the dimension
#'
#' @return A matrix with the encounded restrictions for parameter optimization
#'
#' @examples
#' #p <- 3
#' #compute_ui(p)
compute_ui <- function(p){
  ui   <- matrix(0,nrow = 2 * (p+1) , ncol = (p+1) )
  for(i in 1:(p+1)){
    ui[(2*(i-1) +1),i] <- 1
    ui[(2*(i-1) +2),i] <- -1
  }
  ui
}

#' @rdname compute_ui
compute_ci <- function(p){
  c(-36,-36, rep(c(-10,-10),p ))
}

#' Computes the quantile loss function
#'
#' @param y A double vector with response values
#' @param y_pred A double vector with predictions
#' @param tau A double with values in \eqn{[0,1]} with a probability level to compute the quantile loss function
#'
#' @return Double with the empirical quantile loss
#'
#' @examples
#' quantile_loss(y=2,y=2,.5)
quantile_loss <- function(y, y_pred, tau = .5) {
  error <- y - y_pred
  return(mean((tau - (error < 0)) * error))
}

