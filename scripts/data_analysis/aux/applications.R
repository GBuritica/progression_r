library(engression)
library(grf)

## parameters
prediction_models <- function(X_train,
                              Y_train,
                              X_test,
                              Y_test,
                              nb.iter=10,
                              pred.median=T,
                              fast = T,
                              method = 'llf_progression',
                              q = 0.15,
                              cov = 2){
  models <- list( forest = NULL,
                  engression = NULL,
                  llf = NULL,
                  llf_progression9_mean = NULL,
                  llf_progression9_median = NULL,
                  llf_progression95_mean = NULL,
                  llf_progression95_median = NULL,
                  llf_additive95 = NULL,
                  llf_additive  = NULL
  )
  X_train <- as.matrix(X_train)
  Y_train <- as.vector(Y_train)
  X_test  <- as.matrix(X_test)
  Y_test  <- as.vector(Y_test)
  p  <- ncol(X_train)
  B0 <- 2000
  honesty_use <- TRUE
  ## Fits engression

  ## Fits local linear random forest
  llf <- grf::ll_regression_forest(
    X = X_train,
    Y = Y_train,
    enable.ll.split = TRUE,
    num.trees = B0,
    honesty = honesty_use
  )
  llf_predictions <- predict(llf, X_test )$predictions

  ## Fits random forest
  forest <- grf::regression_forest(
    X = X_train,
    Y = Y_train,
    num.trees = B0,
    honesty = honesty_use
  )
  forest_predictions <- predict(forest, X_test )$predictions

  ## fits progression with a threshold of .9
  if(p==1){
    llf_progression9 <- llf_progression(
      X = X_train,
      Y = Y_train,
      q= q,
      num.trees = B0,
      honesty = honesty_use
    )
    ## mean_progression
    llf_progression9_mean <- predict.llf_progression(llf_progression9, X_test , pred.median = FALSE )$predictions
    ## median_progression
    llf_progression9_median <-predict.llf_progression(llf_progression9, X_test )$predictions
    ## fits progression with a threshold of .95
    llf_progression95 <- llf_progression(
      X = X_train,
      Y = Y_train,
      q= .2,
      num.trees = B0,
      honesty = honesty_use
    )
    ## mean_progression
    llf_progression95_mean  <- predict.llf_progression(llf_progression95, X_test, pred.median = FALSE )$predictions
    ## median_progression
    llf_progression95_median <- predict.llf_progression(llf_progression95, X_test )$predictions
  }
  if(p > 1){
    ##
    ##
    print(cov)
    f_hat <- fit_additive(
      X = X_train,
      Y = Y_train,
      q = q,
      nb.iter = nb.iter,
      plot=FALSE,
      both.sided = TRUE,
      fast = fast,
      cov = cov
    )
    additive_predictions <- f_hat$f(X_test)


    f_hat2 <- fit_additive(
      X = X_train,
      Y = Y_train,
      q = q,
      nb.iter = nb.iter,
      method= 'llf',
      plot=FALSE,
      both.sided = TRUE
    )
    additive_predictions2 <- f_hat2$f(X_test)



  }

  en <- engression::engression(X_train, Y_train, num_epochs = 1000)
  en_predictions <- predict(en, X_test)
  par(mfrow=c(3,2))

  ## results
  models$forest         <- rbind(models$forest, data.frame( 'predictions' = forest_predictions  ) )
  models$engression     <- rbind(models$engression, data.frame( 'predictions' = en_predictions ) )
  models$llf            <- rbind(models$llf, data.frame( 'predictions' = llf_predictions  ) )
  if(p ==1 ){
    models$llf_progression9_mean    <- rbind(models$llf_progression9_mean, data.frame( 'predictions' = llf_progression9_mean) )
    models$llf_progression9_median  <- rbind(models$llf_progression9_median, data.frame( 'predictions' =  llf_progression9_median ) )
    models$llf_progression95_mean   <- rbind(models$llf_progression95_mean, data.frame( 'predictions' = llf_progression95_mean ) )
    models$llf_progression95_median <- rbind(models$llf_progression95_median, data.frame( 'predictions' = llf_progression95_median) )
  }
  if(p > 1){
    models$llf_additive <- rbind(models$llf_additive, data.frame( 'predictions' = additive_predictions) )
    models$llf_additive_llf   <- rbind(models$llf_additive_llf  , data.frame( 'predictions' = additive_predictions2) )
  }


  models
}











