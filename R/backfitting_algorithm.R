#' Backfitting an additive model
#'
#' @description
#' Fits the following additive model using the backfitting algorithm:
#' \deqn{ Y = f_1(X_1) + f_2(X_2) + \cdots f_p(X_p) + \varepsilon_Y}
#' where \eqn{f_1, f_2, \dots f_p} are function such that \eqn{f_i : \mathbb{R} \to \mathbb{R}}
#' and \eqn{\varepsilon_Y} is a noise variable.
#' Each function \eqn{f_i} is fitted separately based on the `method`selected. The default value is `llf_progression` but it also supports random forest and local linear forest fits.
#'
#' @seealso [llf_progression()]
#'
#' @inherit llf_progression
#' @param nb.iter An integer with the number of iterations the back-fitting algorithm should use. The default value equals `1`.
#' @param method A string with values from `rf`, `llf_progression`, `llf` indicating the regression method to use. The default value is `llf_progression`.
#' @param plot Either `TRUE`or `FALSE`to decided whether the fitted regression functions should be plotted.
#' @param reg.fun A list containing the true regression functions \eqn{f_1, f_2, \dots, } used for plotting the regression fit when `plot = TRUE`.
#' @param pred.median A boolean idicating wether the median predictions or mean predictions should be computed
#' @param fast A boolean idicating if the progression method uses the fast implementation
#' @param cov A number indicating the number of covariates to consider. It chooses them automatically based on a measure of extremal depedence.
#' @param ... Additional arguments
#'
#' @return A trained model ready for evaluation.
#' @export
#'
#' @examples
#' set.seed(28)
#' n <- 500
#' reg_fun  <- vector(mode = 'list', length = 2)
#' reg_fun[[1]] <- function(x) (x^2)/3
#' reg_fun[[2]] <- function(x) x
#' X <- matrix(rnorm(n * 2 ), n, 2)
#' f_true <- function(x) reg_fun[[1]](x[, 1]) + reg_fun[[2]](x[, 2])
#' Y <- f_true(X) + rnorm(n)
#' f_hat <- fit_additive(X=X, Y=Y, method = "llf_progression", nb.iter = 5, q = .3 , plot = T)
fit_additive <- function(X, Y,
                         nb.iter = 1,
                         method = 'llf_progression',
                         q = .2,
                         num.trees = 2000,
                         honesty = FALSE,
                         min.node.size = 25,
                         both.sided = TRUE,
                         plot = FALSE,
                         reg.fun = NULL,
                         pred.median=TRUE,
                         fast= TRUE,
                         cov = 2,
                         ...) {
  X <- as.matrix(X)
  p <- ncol(X)
  n <- length(Y)
  a_hat <- stats::median(Y)
  f <- vector("list", p)
  quantile <- .5
  method <- match.arg( method , c('rf', 'llf', 'llf_progression'))
  for (r in 1:nb.iter) { ## iterates backfit
    #print(r)
    ## compute chi coefficient
    if(r == 1){
      xi <- sapply(1:p, function(k) xi_coefficient(X[,k], Y, q = 1-q ) ) ## sorts extremal depedence
      j0 <- sort(-xi, index.return =T)$ix[1]
      p0 <- min(p,cov) ## at most 2 covariates
    }
    #j0 <- 1:p
    ## fits the most extreme signal
    l <- 1
    flag <- TRUE
    while(flag) {
      j <- j0[l]
      #print(j)
      ## computes residuals
      resid_tmp <- Y - a_hat - rowSums(sapply((1:p)[-j], FUN = function(k) {
        if (is.null(f[[k]])) {
          rep(0, times = n)
        } else {
          if (method == "rf"){
            pred    <- utils::getFromNamespace("predict.regression_forest", "grf")
            pred    <- pred(f[[k]], as.matrix(X[, k]) , quantile=quantile )$predictions
            pred  - stats::median(pred)
          }else if (method == 'llf'){
            pred    <- utils::getFromNamespace("predict.ll_regression_forest", "grf")
            if(pred.median){
              pred <- pred( f[[k]], as.matrix(X[, k]) , quantile=quantile )$predictions
              pred - stats::median(pred)
            }
            else{
              pred <- pred(f[[k]], as.matrix(X[, k]) )$predictions
              pred - mean(pred)
            }
          }else if (method == 'llf_progression'){
            pred <- predict.llf_progression(f[[k]], X[, k] , pred.median= pred.median )$predictions
            pred - stats::median(pred)
          }
        }
      }))
      ## looks for next index extremal dependence


      l <- l+1
      if(l>p0){
        flag <- FALSE
      }else{
        if(r==1){
          if(l <= p){
            ind <- which(!(1:p%in%j0))
            xi  <- sapply(ind, function(k) xi_coefficient(X[,k], resid_tmp ) )
            xisort <- sort(-xi, index.return =T)
            if( -xisort$x[1] > 0 ) j0 <- c(j0, ind[sort(-xi, index.return =T)$ix[1]] )
            else flag <- FALSE
            #print(paste0('chi :', -xisort$x))

          }
          else flag <- FALSE
        }

      }
      ## updates residuals data frame

      if (method == "rf") {
        f[[j]] <- grf::regression_forest(
          X = as.matrix(X[, j]),
          Y = resid_tmp,
          num.trees = num.trees,
          honesty = honesty,
          min.node.size = min.node.size
        )

      } else if (method == "llf") {
        f[[j]] <- grf::ll_regression_forest(
          X = as.matrix(X[, j]),
          Y = resid_tmp,
          enable.ll.split = TRUE,
          num.trees = num.trees,
          honesty = honesty,
          min.node.size = min.node.size
        )
      } else if (method == "llf_progression") {
        f[[j]] <- llf_progression(
          X = as.matrix(X[, j]),
          Y = resid_tmp,
          q = q,
          num.trees = num.trees,
          honesty = honesty,
          min.node.size = min.node.size,
          both.sided = both.sided,
          fast = fast
        )
      }
    }

    #print(j0)
    p0 <- length(j0)
    if(plot){
      if(p0 > 1){
        x1_grid <- seq( 2*min(X[,j0[1]]),  2*max(X[,j0[1]]), length.out = 20)
        x2_grid <- seq( 2*min(X[,j0[2]]),  2*max(X[,j0[2]]),length.out = 20)
        ## method for prediction
        if (method == "rf"){
          pred    <- utils::getFromNamespace("predict.regression_forest", "grf")
        }else if (method == 'llf'){
          pred    <- utils::getFromNamespace("predict.ll_regression_forest", "grf")
       }else if (method == 'llf_progression'){
          pred <- predict.llf_progression
        }
       ## graphics
        graphics::par(mfrow=c(1,2))
        if(!is.null(reg.fun)){
          ylim01 <- c(2*min(reg.fun[[j0[1]]](x1_grid)), 2*max(reg.fun[[j0[1]]](x1_grid)))
          ylim02 <- c(2*min(reg.fun[[j0[2]]](x2_grid)), 2*max(reg.fun[[j0[2]]](x2_grid)))
        }
        else{
          ylim01 <- NULL
          ylim02 <- NULL
        }
        graphics::plot(x1_grid, pred(f[[j0[1]]], as.matrix(x1_grid), pred.median = pred.median  )$predictions -
                          stats::median(predict.llf_progression(f[[j0[1]]], as.matrix(X[, j0[1]]), pred.median = pred.median )$predictions),
                     ylab = paste0(
                       'f',j0[1],'(x',j0[1],')') ,
                       xlab = paste0('x',j0[1], '1'),
                     type="l", lty =2,
                     ylim = ylim01,
                     main = paste0('predictions f',j0[1]))
        if(!is.null(reg.fun)) graphics::lines(x1_grid, reg.fun[[j0[1]]](x1_grid))
        graphics::rug(X[,j0[1]], side = 1, col="black")

        graphics::plot(x2_grid, pred(f[[j0[2]]], as.matrix(x2_grid), pred.median =T)$predictions -
                          stats::median(predict.llf_progression(f[[j0[2]]], as.matrix(X[, j0[2]]), pred.median = T )$predictions),
                       ylab = paste0(
                         'f',j0[2],'(x',j0[2],')') ,
                       xlab = paste0('x',j0[2], '1'),
                       type="l", lty =2, col = 'blue',
                       ylim = ylim01,
                       main = paste0('predictions f',j0[2]))
        if(!is.null(reg.fun))  graphics::lines(x2_grid, reg.fun[[j0[2]]](x2_grid), col=
                                               "blue")
        graphics::rug(X[,j0[2]], side = 1, col="blue")
        }
      }

  }
  f_hat <- function(Xnew) {
    a_hat + rowSums(sapply(j0, FUN = function(k) {
      if (is.null(f[[k]])) {
        rep(0, times = n)
      } else {
        if (method == "rf"){
          pred    <- utils::getFromNamespace("predict.regression_forest", "grf")
          pred(f[[k]], as.matrix(Xnew[, k]), quantile=quantile )$predictions - stats::median(pred(f[[k]], as.matrix(X[, k]), quantile=quantile )$predictions)
        }else if (method == 'llf'){
          pred    <- utils::getFromNamespace("predict.ll_regression_forest", "grf")
          if(pred.median) pred(f[[k]], as.matrix(Xnew[, k]), quantile=quantile )$predictions - stats::median(pred(f[[k]], as.matrix(X[, k]), quantile=quantile )$predictions)
          else pred(f[[k]], as.matrix(Xnew[, k]))$predictions - mean(pred(f[[k]], as.matrix(X[, k ]))$predictions)
        }else if (method == 'llf_progression'){
          predict.llf_progression(f[[k]], as.matrix(Xnew[,k]), pred.median = pred.median )$predictions - stats::median(predict.llf_progression(f[[k]], as.matrix(X[, k]), pred.median = pred.median )$predictions)
        }
      }
    }))
  }





  return(f_hat = custom_model(f_hat) )
}


#' Prediction model for additive model
#'
#' Deines the class of the prediction model
#'
#' @param f A list containing the additive functions
#'
#' @return An object of the class "custom_model"
#'
#' @examples
#' #
custom_model <- function(f) {
  model <- list(f = f)
  class(model) <- "custom_model"
  return(model)
}

#' Prediction function for additive model
#'
#' Defines the function for predicting an output from the custom_model
#'
#' @param object An object of class "custom_model".
#' @param newdata A matrix with \eqn{p} columns, where \eqn{p} is the number of predictors from the model.
#' @param ... Additional arguments.
#'
#' @return The evaluation of the additive model
#'
#' @examples
#' #
predict.custom_model <- function(object, newdata, ...) {
  result <- object$f(newdata)
  return(result)
}






#' Compute coefficients of extremal dependence
#'
#' Method for calculation the xi coefficient of extremal dependence.
#'
#' @param data1 A vector of length \eqn{n}.
#' @param data2 A vector of same length as `data1`.
#' @param q A probability of exceedances.
#'
#' @return The value of the xi coefficient.
#'
#' @examples
#' #
xi_coefficient <- function(data1,data2, q = .9 ){
  data1 <- abs(data1)
  data2 <- abs(data2)
  ind <-  ((data2 > stats::quantile(data2, q)) | (data2 < stats::quantile(data2, 1-q)) ) & ((data1 > stats::quantile(data1, q)) | (data1 < stats::quantile(data1, 1-q)) )
  mean(ind)/(2*(1-q))
}

