### Check on several examples the fit for each model
figure_dest_folder <- 'results/plots'
data_dest_folder <- 'results/data'
library(progression)
library(engression)
library(grf)
library(ggplot2)
source('scripts/plots/save_plot.R')
## Univariate example


#   Setting 1
set.seed(123)
reg_fun <- function(x) x^3 / 10
df_X <- 3
k <- 2
train_dist <- function(n) k*rt(n, df = df_X)
eps_dist <- function(n) 4 * rnorm(n)


# # Setting 2
df_X <- 3
k <- 4
reg_fun <- function(x) ifelse(x > 0, (1+x)^{1/2} -1  , -(1-x)^{1/2} + 1)
train_dist <- function(n) k*rt(n, df = df_X)
eps_dist <-  function(n) .5*rnorm(n)



# Setting 3
set.seed(28)
df_X <- 3
k <- 1
reg_fun <- function(x) x^2/8
train_dist <- function(n) k*rt(n, df = df_X)
eps_dist <-  function(n) .1*rnorm(n)


engress = T # Set engression to false to for a faster implementation


n_rep <- 5
n <- 1000
n_test <- 50
qu <- k*qt(1 - 1 / (4 * n), df = df_X) # max(abs(X_train)) * 1.3 #60
X_test <- seq(-qu, qu, length = n_test)
Y_true <- reg_fun(X_test)

q <- 0.3
honesty_use <- FALSE
B0 <- 2000

# Set the dimensions of boxplot data
if(engress) dim_data <- c(5, n_rep, n_test) # 5 methods, n_rep repetitions, n_test predictions
if(!engress) dim_data <- c(4, n_rep, n_test) # 4 methods, n_rep repetitions, n_test predictions
res_sim <- array(NA, dim = dim_data)


for (r in 1:n_rep) {
  print(r)

  X_train <- train_dist(n)
  eps_Y <- eps_dist(n)
  Y_train <- reg_fun(X_train) + eps_Y

  # Fit engression
  if(engress) en <- engression(X_train, Y_train, num_epochs = 1000)

  ## Fit at locally linear random forest
  llf <- ll_regression_forest(
    X = as.matrix(X_train),
    Y = Y_train,
    enable.ll.split = TRUE,
    num.trees = B0,
    honesty = honesty_use
  )


  ## Fit at random forest
  forest <- regression_forest(
    X = as.matrix(X_train),
    Y = Y_train,
    num.trees = B0,
    honesty = honesty_use
  )


  ## fits fast progression
  llf_LP  <- llf_progression(
    X = as.matrix(X_train),
    Y = Y_train,
    q = q,
    honesty = honesty_use,
    both.sided = TRUE,
    fast = TRUE
  )
  pred <- predict(llf_LP, as.matrix(X_test), pred_median = TRUE)$predictions


  ## fits slow progression
  llf_LP_both_sided <- llf_progression(
    X = as.matrix(X_train),
    Y = Y_train,
    q = q,
    honesty = honesty_use,
    both.sided = TRUE,
    fast = FALSE
  )


  llf_pred <- predict(llf, as.matrix(X_test))$predictions
  forest_pred <- predict(forest, as.matrix(X_test))$predictions
  llf_prog_median_pred <- pred
  llf_prog_median_both_pred <- predict(llf_LP_both_sided, as.matrix(X_test), pred_median = TRUE, both_sided = TRUE)$predictions
  if(engress) {
    en_pred <- predict(en, X_test)
    res_sim[, r, ] <- rbind(
      forest_pred, llf_pred, llf_prog_median_pred,
      llf_prog_median_both_pred, en_pred
    )
  }
  if(!engress) {
    res_sim[, r, ] <- rbind(
      forest_pred, llf_pred, llf_prog_median_pred,
      llf_prog_median_both_pred
    )
  }


}

#save(res_sim, file = 'cubic.Rdata' )
#load("figures/res_sim_sqrt.Rdata")
library(matrixStats)
alpha <- .2
#pdf("sqrt_all.pdf", width = 10, height = 3)
if(engress == F) par(mfrow=c(1,4))
if(engress == T) par(mfrow=c(1,5))
method_idx <- 1
colors_with_alpha <- rgb(0, 0, 0, alpha = alpha)
matplot(X_test, t(rbind(Y_true, res_sim[method_idx, , ])), type = "l", lty = 1, col = colors_with_alpha, lwd = 1, ylim=1*range(Y_true), main="Random forest", ylab = '', xlab = 'x test')
points(X_train, Y_train, col = colors_with_alpha, cex = 0.2)
lines(X_test, Y_true, col = "blue", lwd = 2, lty=1)
lines(X_test, rowMedians(t(rbind(Y_true, res_sim[method_idx, , ]))), col="orange", lwd=2)

method_idx <- 2
colors_with_alpha <- rgb(0, 0, 0, alpha = alpha)
matplot(X_test, t(rbind(Y_true, res_sim[method_idx, , ])), type = "l", lty = 1, col = colors_with_alpha, lwd = 1, ylim=1*range(Y_true), main="Local linear forest", ylab = '', xlab = 'x test')
points(X_train, Y_train, col = colors_with_alpha, cex = 0.2)
lines(X_test, Y_true, col = "blue", lwd = 2)
lines(X_test, rowMedians(t(rbind(Y_true, res_sim[method_idx, , ]))), col="orange", lwd=2, ylab = '', xlab = 'x test')

method_idx <- 3
colors_with_alpha <- rgb(0, 0, 0, alpha = alpha)
matplot(X_test, t(rbind(Y_true, res_sim[method_idx, , ])), type = "l", lty = 1, col = colors_with_alpha, lwd = 1, ylim=1*range(Y_true), main="Progression", ylab = '', xlab = 'x test')
points(X_train, Y_train, col = colors_with_alpha, cex = 0.2)
lines(X_test, Y_true, col = "blue", lwd = 2)
lines(X_test, rowMedians(t(rbind(Y_true, res_sim[method_idx, , ]))), col="orange", lwd=2, ylab = '', xlab = 'x test')

if(engress == T){
  method_idx <- 5
  colors_with_alpha <- rgb(0, 0, 0, alpha = alpha)
  matplot(X_test, t(rbind(Y_true, res_sim[method_idx, , ])), type = "l", lty = 1, col = colors_with_alpha, lwd = 1, ylim=1*range(Y_true), main="Engression", ylab = '', xlab = 'x test')
  points(X_train, Y_train, col = colors_with_alpha, cex = 0.2)
  lines(X_test, Y_true, col = "blue", lwd = 2)
  lines(X_test, rowMedians(t(rbind(Y_true, res_sim[method_idx, , ]))), col="orange", lwd=2, ylab = '', xlab = 'x test')

}

#dev.off()

method_idx <- 4
colors_with_alpha <- rgb(0, 0, 0, alpha = alpha)
matplot(X_test, t(rbind(Y_true, res_sim[method_idx, , ])), type = "l", lty = 1, col = colors_with_alpha, lwd = 1, ylim=1*range(Y_true), main="cautiousProgression", ylab = '', xlab = 'x test')
points(X_train, Y_train, col = colors_with_alpha, cex = 0.2)
lines(X_test, Y_true, col = "blue", lwd = 2)
lines(X_test, rowMeans(t(rbind(Y_true, res_sim[method_idx, , ]))), col="orange", lwd=2, ylab = '', xlab = 'x test')

#sum(apply(res_sim[2, , ],  1, function(x)  mean(abs( x- Y_true)) ))
#sum(apply(res_sim[4, , ],  1, function(x)  mean(abs( x- Y_true)) ))






