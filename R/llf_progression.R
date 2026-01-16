#' Local linear progression forest
#'
#' @description Trains local linear progression for predicting a response variable \eqn{Y} based on covariates \eqn{X}.
#' It estimates by default a conditional median, but can also be used for quantile regression.
#' The method tackles accurate prediction on new test points \eqn{x_{\text{test}}} which are extreme or underrepresented during training.
#' To do so, it carries out a standarization of margins to the Laplace distribution scale.
#' Hence, no previous marginal standarization is needed.
#'
#' This function calls the [grf::ll_regression_forest()] and some of its parameters are inherited from this function.
#'
#' @seealso [grf::ll_regression_forest()]
#'
#' @param X A double vector, matrix or data.frame containing covariates.
#' @param Y A vector of responses with continuous outcomes.
#' @param q A probability taking values in \eqn{(0,1/2)}, indicating the probability of threshold exceedances. In this case, \eqn{q} percentage of the data is used for the GPD fit during marginal transformation. The default value is \eqn{q = 0.3}.
#' @param num.trees An integer indicating the number of trees in the forest. Default is 2000.
#' @param honesty A boolean of whether to use honest splitting or not. Default is FALSE.
#' @param both.sided A boolean. Default value is `FALSE`, and currently this feature is only supported for a univariate predictor \eqn{X}. If `TRUE` then the marginal standarizations of \eqn{X} and \eqn{Y} are computed separately for upper data: \eqn{\{(X,Y) | X > median(X)\}} and lower data: \eqn{\{(X,Y) | X \leq median(X)\}}. The feature `both.sided = TRUE` is convenient when the upper and lower tails of \eqn{Y | X > median(X) } and \eqn{Y | X \leq median(X) } are not equivalent.
#' @param min.node.size An integer with the target for the minimum number of observation in each tree leaf. Default is 5.
#' @param fast A boolean that when is set to false returns a conservative extrapolation. It computes a conservative response based on profile likelihood confidence intervals for the generalized Pareto distribution fit
#' @param FX A function with the cumulative distribution function (CDF) of the predictor X.
#' @param FY A function with the CDF of the response Y.
#' @param QX A function with the quantile function of the predictor X.
#' @param QY A function with the quantile function of the response Y.
#'
#' @return A trained local linear progression object
#' @export
#'
#' @examples
#' # Train a standard progression object
#' n <- 50
#' p <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] + rnorm(n)
#' llf_progression_forest <- llf_progression(X, Y)
llf_progression <- function(X,
                            Y,
                            q = 0.1,
                            num.trees = 2000,
                            honesty = FALSE,
                            both.sided = TRUE,
                            min.node.size = 25,
                            fast = FALSE,
                            FX = NULL,
                            FY = NULL,
                            QX = NULL,
                            QY = NULL
) {

  n <- length(Y)
  X <- as.matrix(X)
  p <- ncol(X)
  # lower and upper threshold
  upper_X = NULL
  lower_X = NULL
  upper = NULL
  # indices for large and small X values when p = 1
  up = NULL
  low = NULL
  mid = NULL
  # probability level
  q <- min( q, n^(0.7-1)) ## probability level
  if(!is.null(FX)) both.sided = FALSE
  if(p > 1){
    if(both.sided) warning('using both.sided = FALSE, as currently this option is only supported for a univariate covariate')
    both.sided = FALSE
  }else{
    if(both.sided){

      upper_X = stats::quantile(X[,1],(1-q)) ## used to separate data above and below threshold
      lower_X = stats::quantile(X[,1],q)
      up = (X[,1] > upper_X)
      low = (X[,1] < lower_X)
      mid   = !(up + low)
      upper = list(upper=up, lower = low, mid = mid)
    }
  }

  ## Laplace transformation applied separately on different low, bulk and high X values.
  trans_X <- apply(X,2,function(k) prepare_laplace_transformations(data = k, upper=upper, fast = fast, cdf = FX, Qf = QX ))
  trans_Y <- prepare_laplace_transformations(data = Y,  upper=upper, fast = fast, cdf = FY, Qf = QY )

  X_train_LP <- NULL
  for(i in 1:p) X_train_LP <- cbind(X_train_LP, trans_X[[i]]$data_to_laplace(X[,i], upper=upper))
  Y_train_LP <- trans_Y$data_to_laplace(Y,upper=upper)

  if(both.sided){
    llf_LP <- list()
    llf_LP$upper <- grf::ll_regression_forest(
      X = as.matrix(X_train_LP[up]),
      Y = Y_train_LP[up],
      enable.ll.split = TRUE,
      num.trees = num.trees,
      honesty = honesty,
      min.node.size = min.node.size
    )
    llf_LP$lower <- grf::ll_regression_forest(
      X = as.matrix(X_train_LP[low]),
      Y = Y_train_LP[low],
      enable.ll.split = TRUE,
      num.trees = num.trees,
      honesty = honesty,
      min.node.size = min.node.size
    )
    llf_LP$mid <- grf::ll_regression_forest(
      X = as.matrix(X_train_LP[mid]),
      Y = Y_train_LP[mid],
      enable.ll.split = TRUE,
      num.trees = num.trees,
      honesty = honesty,
      min.node.size = min.node.size
    )
    graphics::par(mfrow=c(1,1))
    X_train_LP <- list( upper = as.matrix(X_train_LP[up,]), lower =  as.matrix(X_train_LP[low,]), mid = as.matrix(X_train_LP[mid,]) )
    Y_train_LP <- list( upper = Y_train_LP[up], lower =  Y_train_LP[low], mid = Y_train_LP[mid])
    #plot(Y_train_LP$upper, predict(llf_LP$upper, as.matrix(X_train_LP$upper))$predictions )
    #abline(0,1)
  }else{
    llf_LP <- grf::ll_regression_forest(
      X = as.matrix(X_train_LP),
      Y = Y_train_LP,
      enable.ll.split = TRUE,
      num.trees = num.trees,
      honesty = honesty,
      min.node.size = min.node.size
    )
  }

  model <- list(llf_LP = llf_LP, trans_X = trans_X, trans_Y = trans_Y, X_train_LP = X_train_LP, Y_train_LP = Y_train_LP, both.sided = both.sided, upper_X = upper_X, lower_X = lower_X, upper = upper)
  class(model) <- "llf_progression"
  return(model)

}


