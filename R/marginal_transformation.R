#' computes semi-parametric CDF and quantile functions
#'
#' This function return both the CDF and the quantile function of the data based on a semi-parametric approach. In the bulk of the distribution, the method used a non-parametric approximation, while in the tails it used the Generalized Pareto distribution (GPD) fit
#'
#' @param data A double vector corresponding to the data from which the semi-parametric CDF will be computed
#' @param fast A binary variable indicating if a profile likelihood grid search or a fast point estimate must be returned for the quantile
#' @param method A String containing either 'MLE' or 'proflik' to define how confidence intervals are computed for the conservative extrapolation implementation.
#' @param cdf A function containing the cumulative distribution function of the data or 'NULL' if this should be computed by the method.
#' @param Qf A function with the quantile function of the data or 'NULL' if this function should be computed by the method.
#' @param p A small probability indicating the proportion of data considered as extremes in the tails.
#'
#' @return A list containing two functions:
#'    * `cdf_func(x)` : The CDF function of the data taking values from \eqn{\mathbb{R}} to \eqn{(0,1)}.
#'    * `quantile_func(q)` : The quantile function of the data taking values from \eqn{(0,1)} to \eqn{\mathbb{R}}.
#' @export
#'
#' @examples
#' set.seed(28)
#' data <- rnorm(100)  # Generate some normal data
#' functions <- prepare_semi_parametric_functions(data, method='proflik', fast=T, p = 0.2)
#' curve(functions$cdf(x), from = -4, to = 4, main = 'CDF')
#' p <- functions$cdf(data)
#' q <- functions$inverse_cdf(p)
#' plot(data,q)
#' abline(0,1)
prepare_semi_parametric_functions <-
  function(data,
           fast = TRUE,
           method = 'MLE',
           cdf = NULL,
           Qf = NULL,
           p = 0.5) { ## fits upper and lower tails with GPD
    lambda <- 0
    n <- length(data)
    sorted_data <- sort(data)

    lower_threshold <-
      stats::quantile(data, p)          # compute lower quantile
    upper_threshold <- stats::quantile(data, 1 - p)
    emp <- stats::ecdf(data)
    # Fit GPD to the tails
    upper_tail <- data[data > upper_threshold]
    lower_tail <- data[data < lower_threshold]
    upper_fit <- mev::fit.gpd(upper_tail - upper_threshold, 0)
    lower_fit <- mev::fit.gpd(lower_threshold - lower_tail, 0)

    # end points
    lower_end_point <-
      ifelse(lower_fit$param[2] < 0,
             lower_threshold + lower_fit$param[1] / lower_fit$param[2],
             -Inf) ## modif. end points
    upper_end_point <-
      ifelse(upper_fit$param[2] < 0,
             upper_threshold - upper_fit$param[1] / upper_fit$param[2],
             Inf) ## modif. end points
    # CDF Function
    if(is.null(cdf)){
      cdf_func <- function(x) {
        if (x < lower_threshold) {
          # Lower tail GPD CDF
          est <-
            p * (
              1 - evd::pgpd(
                lower_threshold - x,
                loc = 0,
                scale = lower_fit$param[1],
                lower_fit$param[2]
              )
            )
          return(est)
        } else if (x > upper_threshold) {
          # Upper tail GPD CDF
          est <-
            1 - p * (
              1 - evd::pgpd(
                x - upper_threshold,
                loc = 0,
                scale = upper_fit$param[1],
                shape = upper_fit$param[2]
              )
            )
          return(est)
        } else {
          # Empirical CDF for the central part
          return(emp(x))
        }
      }
    }else cdf_func <- function(x) cdf(x)

    # Inverse CDF Function
    if(is.null(Qf)){
      quantile_func <- function(q) {
        if (q < p) {
          # Lower tail inverse GPD
          q0   <-
            sapply(round(
              1 - round(round(q, 14)/ p,14) , 14) , function(k)
                min(k, 0.9994))#
          res2 <-
            function(q0)
              lower_threshold - evd::qgpd(
                q0,
                loc = 0,
                scale = lower_fit$param[1],
                shape = lower_fit$param[2]
              )
          if (any(is.na(lower_fit$std.err)))
            fast = TRUE
          if (!fast) {
            se_lower <- stats::confint(lower_fit)
            if (method == 'MLE') {
              ### normal ci
              fit <-
                extRemes::fevd(
                  lower_threshold - lower_tail ,
                  threshold = 0,
                  type = 'GP' ,
                  time.units = 'years' ,
                  method = 'MLE'
                )
              res <- function(q0) {
                res0 <- extRemes::ci.fevd.mle(
                  fit,
                  type = 'return.level',
                  return.period = 1 / (1 - q0),
                  alpha = 0.8
                )    ## It is  possible to do proflik here, but convergence problems might arise.
                #data.frame(p = q0, estimate = unname(res0[2]), low_ci = unname(res0[1]), up_ci = unname(res0[3]) )
                lower_threshold - unname(res0[1])
              }

            }
            if (method == 'proflik') {
              res <- function(q0) {
                xi <- mev::gpd.pll(
                  psi = seq(
                    from = evd::qgpd(
                      q0,
                      loc = 0,
                      scale = se_lower[1, 1],
                      shape = se_lower[2, 1]
                    ),
                    to = evd::qgpd(
                      q0,
                      loc = 0,
                      scale = se_lower[1, 2],
                      shape = se_lower[2, 2]
                    ) ,
                    length.out = 20
                  ),
                  param = 'quant',
                  p = 1 - q0,
                  mod = "profile",
                  dat = lower_threshold - lower_tail,
                  correction = TRUE,
                  threshold = 0,
                  plot = FALSE
                )
                con_xi <-
                  mev:::confint.eprof(xi,
                                      level = 0.5,
                                      method = 'smooth',
                                      spar = 0.7)
                lower_threshold - unname(con_xi[2])
              }
            }
          }
          if (fast)
            res <- res2

          return(sapply(q0, res))
        } else if (q > 1 - p) {
          # Upper tail inverse GPD
          q0   <-
            sapply((1 - round(
              round(1 - q, 14) / p, 14)), function(k)
              min(k, 0.9994))
          #res <- function(q0) upper_threshold + evd::qgpd(q0, loc = 0, scale = upper_fit$param[1], shape = upper_fit$param[2])

          res2 <-
            function(q0)
              upper_threshold + evd::qgpd(
                q0,
                loc = 0,
                scale = upper_fit$param[1],
                shape = upper_fit$param[2]
              )
          if (any(is.na(upper_fit$std.err)))
            fast = T
          if (fast)
            res <- res2
          else{
            se_upper <- stats::confint(upper_fit)
            if (method == 'MLE') {
              ### normal ci
              fit <-
                extRemes::fevd(
                  upper_tail - upper_threshold ,
                  threshold = 0,
                  type = 'GP' ,
                  time.units = 'years' ,
                  method = 'MLE'
                )
              res <- function(q0) {
                res0 <- extRemes::ci.fevd.mle(
                  fit,
                  type = 'return.level',
                  return.period = 1 / (1 - q0),
                  alpha = 0.8
                )    ## It is  possible to do proflik here, but convergence problems might arise.
                #data.frame(p = q0, estimate = unname(res0[2]), low_ci = unname(res0[1]), up_ci = unname(res0[3]) )
                upper_threshold + unname(res0[1])
              }


            }
            if (method == 'proflik') {
              res <- function(q0) {
                xi <- mev::gpd.pll(
                  psi = seq(
                    from = evd::qgpd(
                      q0,
                      loc = 0,
                      scale = se_upper[1, 1],
                      shape = se_upper[2, 1]
                    ),
                    to = evd::qgpd(
                      q0,
                      loc = 0,
                      scale = se_upper[1, 2],
                      shape = se_upper[2, 2]
                    ) ,
                    length.out = 20
                  ),
                  param = 'quant',
                  p = 1 - q0,
                  mod = "profile",
                  dat = upper_tail - upper_threshold,
                  correction = TRUE,
                  threshold = 0,
                  plot = FALSE
                )
                con_xi <-
                  mev:::confint.eprof(xi,
                                      level = 0.5,
                                      method = 'smooth',
                                      spar = 0.7)
                upper_threshold + unname(con_xi[2])
              }
            }

          }
          return(sapply(q0, res))
        } else {
          # Empirical quantile for the central part
          central_quantile_index <- floor(q * (n + 1))
          return(sorted_data[central_quantile_index])
        }
      }
    }else quantile_func <- function(q) Qf(q)

    return(list(
      cdf = Vectorize(cdf_func),
      inverse_cdf = Vectorize(quantile_func)
    ))
  }



#' computes Laplace and inverse Laplace transformations
#'
#' @description
#' This function returns estimates of the Laplace and inverse Laplace transformations based on the data. Remember that the Laplace transformation of a random variable \eqn{X} is given by
#' \deqn{L_X  = Q_L \circ F_X}
#' where \eqn{F_X} is the CDF of \eqn{X} and \eqn{Q_L} is the quantile function of the Laplace distribution.
#'
#' @seealso [prepare_semi_parametric_functions()] for estimation of the CDF function \eqn{F_X}.
#' @inheritParams prepare_semi_parametric_functions
#' @param upper Either NULL or a boolean vector with the length of the data. The default value is NULL meaning the Laplace transformation is with respect to all the data.  Otherwise, the function computes separately the CDF on the two sets: `data[upper]` and `data[!upper]`.
#'
#' @return A list containing two functions
#'    * `data_to_laplace(x, upper = NULL)` : a mapping taking a double vector x to its Laplace-transformed values \eqn{L_X(x)}. In the case the upper and lower data Laplace transformations were computed, then it also takes `upper` as an input. It should be a boolean vector with the same length as $x$ indicating whether the values of \eqn{x} belong to the upper side of the data or not.
#'    * `laplace_to_data(x, upper = NULL)` : mapping taking a double vector x to its inverse Laplace-transformed values. Similarly as for the `data_to_laplace` function, it can take as an argument the boolean vector `upper`. This is a marker of the \eqn{x} values belonging to the upper of lower part of the data.
#' @export
#'
#' @examples
#' set.seed(28)
#' data <- rt(1000, df = 2)
#' transformations <- prepare_laplace_transformations(data, p = 0.2)
#' laplace_values <- transformations$data_to_laplace(data)
#' reconstructed_data <- transformations$laplace_to_data(laplace_values)
#' plot(data, laplace_values, main = "Data vs. Laplace Transformation")
#' plot(reconstructed_data, laplace_values,  main = "Data Reconstruction vs. Laplace Transformation")

prepare_laplace_transformations <- function(data, upper = NULL,  fast =TRUE, cdf = NULL, Qf = NULL, p = 0.2 ) { ## upper is the set of indices where the cdf should be computed
  # Fit the GPD and prepare the CDF and inverse CDF functions
  if (!is.null(upper)) {
    functions_upper <- prepare_semi_parametric_functions(data[upper$upper], fast=fast) ## fits semi parametric to data above the 1-pth quantile
    functions_lower <- prepare_semi_parametric_functions(data[upper$lower],  fast=fast) ## fits semi parametric to data below the pth quantile
    functions_mid   <- prepare_semi_parametric_functions(data[upper$mid],   fast=fast) ## fits semi parametric to data in between the upper and lower quantiles
  } else {
    functions_upper <- prepare_semi_parametric_functions(data,  fast, cdf = cdf, Qf = Qf, p = p) ## fits semi-parametric to all data
    functions_lower <- functions_upper
    functions_mid   <- functions_upper
  }

  # Function to transform data scale to Laplace scale
  data_to_laplace <- function(x, upper = NULL) {
    if (is.null(upper)) upper <- list( upper = rep(TRUE, length(x)), lower = rep(TRUE, length(x)), mid = NULL)
    # Semi-parametric CDF to get probabilities
    prob <- ifelse(upper$upper,functions_upper$cdf(x), functions_mid$cdf(x))
    prob <- ifelse(upper$lower,functions_lower$cdf(x), prob )
    # Transforming probability to Laplace scale assuming mu = 0, b = 1
    laplace_values <- ifelse(prob >= 0.5,
                             -log(2 * (1 - prob)),
                             log(2 * prob)
    )
    return(unname(bound_exp(laplace_values)))
  }

  # Function to transform Laplace scale to data scale
  laplace_to_data <- function(y, upper = NULL) {
    # Convert Laplace values to probabilities assuming mu = 0, b = 1
    if (is.null(upper)) upper <- list( upper = rep(TRUE, length(y)), lower = rep(TRUE, length(y)), mid = NULL)
    prob <- ifelse(y >= 0,
                   1 - 0.5 * exp(-y),
                   0.5 * exp(y)
    )
    # Semi-parametric inverse CDF to get data values
    data_values <- ifelse(upper$upper, functions_upper$inverse_cdf(prob), functions_mid$inverse_cdf(prob))
    data_values <- ifelse(upper$lower, functions_lower$inverse_cdf(prob), data_values)
    return(unname(data_values))
  }

  return(list(
    data_to_laplace = data_to_laplace,
    laplace_to_data = laplace_to_data
  ))
}
