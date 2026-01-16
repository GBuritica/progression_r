#' Spline generator from a B-basis
#'
#' genSpline is a generator of a polynomial splines based on the B-slpline basis matrix
#'
#' @param x A double vector with the values to evaluate
#' @param knots A double vector with breakpoints defining the spline.
#' @param degree An integer with the polynomial degree.
#' @param theta A double vector of length equal to df = length(knots) + degree + 1
#'
#' @return A data.table with the spline attributes
#'
genSpline <- function(x, knots, degree, theta) {
  basis <- splines::bs(x = x, knots = knots, degree = degree,
                       Boundary.knots = c(0,1), intercept = TRUE)

  y.spline <- basis %*% theta

  dt <- data.table::data.table(x, y.spline = as.vector(y.spline))
  return(list(dt = dt, basis = basis, knots = knots))
}


#' Spline function
#'
#' @description Evaluates a vector on the spline function
#'  * `spline1()` returns the values of one spline
#'  * `spline2()` returns the values of a different spline
#'  * `spline3()` different spline
#'  * `spline4()` different spline
#'  * `spline5()` different spline
#'  * `spline6()` different spline
#'  * `spline7()` different spline
#'  * `spline8()` different spline
#'
#'
#' @param x A double vector where the spline will be evaluated
#' @param qua An integer for scaling the y-axis output. The default value is qua = 20.
#' @param degree A boolean indicating weather the output should be a list with the spline values together with the function degree/exponent. The default value is `FALSE`.
#'
#' @returns
#'
#'  * If degree equals `FALSE` : a double vector of same length as x with the spline values.
#'  * If degree equals `TRUE` : a list of two elements : the spline values, and a double value with the spline degree or function exponent
#' @export
#'
#' @examples
#' library(progression)
#' x <- seq(-30,30,by = .1)
#' y <- spline1(x)
#' plot(x,y, type = 'l')
spline1 <- function(x, qua = 15, degree=F){
  knots <- c(0.3, 0.5, 0.7)
  theta <- c(0.6, 0.1, 0.5, 0.2, 0.8, 0.3)
  gammaf <- 2
  y <-    2  * genSpline( .5 + x/qua, knots, gammaf, theta)$dt$y.spline
  if(degree) return(list(y=y,degree=gammaf))
  else y
}



#' @rdname spline1
spline2 <- function(x, qua = 15, degree=F){
  knots <- c(0.4, 0.6)
  theta <- c(0.2, 0.4, 0.1, 0.9, 0.6)
  gammaf <- 2
  y <-  4 * genSpline(.5 + x/qua, knots, gammaf, theta)$dt$y.spline
  if(degree) return(list(y=y,degree=gammaf))
  else y
}


#' @rdname spline1
spline3 <- function(x, qua = 15, degree=F){
  knots <- c(0.4, 0.6)
  theta <- c(0.2, 0.6, 0.1, 0.9, 0.2, 0.8)
  gammaf <- 3
  y <- 1 * genSpline( .5 +x/qua, knots, gammaf, theta)$dt$y.spline
  if(degree) return(list(y=y,degree=gammaf))
  else y
}


#' @rdname spline1
spline4 <-  function(x, qua = 15, degree = F){
  y <-   2 * sign(x)*abs(2*x)^.2
  gammaf <- .2
  if(degree) return(list(y=y,degree=gammaf))
  else y
}

#' @rdname spline1
spline5 <- function(x, qua = 15, degree = F){
  y <- - (1/20) * x + (1/500) * x^3
  #y <- (1/5) * x^3
  gammaf <- 3
  if(degree) return(list(y=y,degree=gammaf))
  else y
}


#' @rdname spline1
spline6  <- function(x, qua = 15, degree = F){
  y <-   2 * sign(x)*abs(2*x)^.5
  gammaf <- .5
  if(degree) return(list(y=y,degree=gammaf))
  else y
}

#' @rdname spline1
spline7  <- function(x, qua = 15, degree = F){
  y <-   8 * abs(x/qua)^2
  gammaf <- 2
  if(degree) return(list(y=y,degree=gammaf))
  else y
}

#' @rdname spline1
spline8  <- function(x, qua = 15, degree = F){
  y <-  10 * sign(x)*abs(x/qua)^4
  gammaf <- 4
  if(degree) return(list(y=y,degree=gammaf))
  else y
}


#' @rdname spline1
spline9  <- function(x, qua = 15, degree = F){
  y <- 8 * log(20 + abs(x))
  gammaf <- 1
  if(degree) return(list(y=y,degree=gammaf))
  else y
}




