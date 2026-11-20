## We follow Brockwell and Davis to remove Trend and Season
## Remove Trend
## k is the number of day lags, k > 1
## d is the number of years
## data is a data_frame whose first column is the date
adjust_mean         <- function(x,k) sapply((k+1):length(x), function(i) mean(x[(i-k):(i-1)], na.rm = T))
## Season 1 year
## d number of years
adjust_season       <- function(x,d,year=365.25)  sapply(  (floor(d*year)+1):length(x),  ## fijo el dia
                                                           function(i) seasonal_mean(  sapply(1:d, function(k) concatenate_years(x,i,k,year) ) ) ) ## valor de los últimos 4 años
concatenate_years   <- function(x,i,k,year) x[ (i- floor(k*year) ) ]

seasonal_mean       <- function(x){
  if(!all(is.na(x)) ) mean( x, na.rm=T )
  else NA
}

stationarity_adjust <- function(data, d = 10,k = 20,years = 365.25, date=TRUE){

  ## adjust trend
  if(date) data[,1] <- as.Date(data[,1])
  data_mean_adjust       <- data
  n                      <- length(data[1,]) ## number of columns
  data_trend             <- apply(data[,2:n], 2, function(x) adjust_mean(x,k)  )
  data_mean_adjust[-(1:k),2:n]  <- data_mean_adjust[-(1:k),2:n] - data_trend
  data_trend                    <- as.data.frame(data_trend)
  if(date) data_trend$Date               <- as.Date(data_mean_adjust[-(1:k),1])
  else     data_trend$Date               <- (data_mean_adjust[-(1:k),1])

  ## adjust season
  data_season    <- apply(data_mean_adjust[,2:n], 2, function(x) adjust_season(x,d,years)  )

  data_adjusted <- data
  data_adjusted[-(1:floor(d*years)),2:n ] <- data_mean_adjust[-(1:floor(d*years)),2:n ]  - data_season
  data_season      <- as.data.frame(data_season)
  if(date) data_season$Date <- as.Date(data_adjusted[-(1:floor(d*years)),1] )
  else     data_adjusted[-(1:floor(d*years)),1]


  adjusted <- list( 'data_trend' = data_trend,
                    'data_season' = data_season,
                    'data_adjusted' = data_adjusted,
                    'data_original' = data,
                    'd' = d,
                    'years' = years,
                    'k' = k )
  return(adjusted)
}

## data should start with a date
block_maxima <- function(data, b = 8, fun = max ){

  n     <- length(data[,1])
  p     <- length(data[1,])
  names <- colnames(data)
  m     <- floor(n/b)
  date  <- as.Date(sapply(1:m, function(k) data[ ((k-1)*b +1), 1 ] ))
  data2 <- sapply(2:p, function(l) sapply( 1:m, function(k) fun(data[( ((k-1)*b+1):(k*b) ) , l] ) ) )
  data2 <- cbind(date,data2)
  data2 <- as.data.frame(data2)
  colnames(data2) <- names
  return(data2)

}

## for example
#precip_adjusted <- stationarity_adjust(precip)


###########################################################
## Adding noise to remove zeros
## It must receives a time series
##
add_noise_marginal  <- function(x) return(x + rnorm(n=length(x), sd=stats::sd(x,na.rm=T)/2))
add_noise <- function(data){
  data[,2:length(data[1,])] <- apply(data[,2:length(data[1,])], 2, add_noise_marginal  )
  return(data)
}

