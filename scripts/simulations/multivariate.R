source("scripts/simulations/functions.R")
# devtools::install_local("../progression", force=TRUE)

set.seed(12144)


# sim_setting <- "polynomial" #"Friedman"

for(sim_setting in c("polynomial", "Friedman")){


  if (sim_setting == "polynomial") {
    reg_fun1 <- function(x) {ifelse(x > 0, 4 * abs(x)^{1 / 2}, -4 * abs(x)^{1 / 2})}
    reg_fun2 <- function(x) x^3 / 20
    rho <- 0 ## can adjust the dependence between predictors
    train_dist <- function(n) qt(apply(rmvnorm(n, sigma = cbind(c(1, rho), c(rho, 1))), 2, unif), df = 4)
    f_true <- function(x) reg_fun1(x[, 1]) + reg_fun2(x[, 2])
    eps_dist <- function(n) 2 * rnorm(n)
  }


  # # Setting 3: Friedman function
  if (sim_setting == "Friedman") {
    train_dist <- function(n) 1 / 5 * qt(apply(rmvnorm(n, mean = rep(0, 5)), 2, unif), df = 4)
    reg_fun1 <- function(x) 20 * (x - 0.5)^2
    reg_fun2 <- function(x) 10 * x
    reg_fun3 <- function(x) 5 * x
    f_true <- function(x) reg_fun1(x[, 1]) + reg_fun2(x[, 2]) + reg_fun3(x[, 3]) + 10 * sin(pi * x[, 4] * x[, 5])
    eps_dist <- function(n) 2 * rnorm(n)
  }


  n_rep <- 2

  n <- 1000

  n_test <- 200
  q <- 0.1
  nb.iter <- 2
  honesty_use <- FALSE
  B0 <- 2000

  # Set the dimensions of boxplot data
  dim_data <- c(4, 4, n_rep, 3)  # 4 methods, 4 shifts, n_rep repetitions, 3 error measures
  res_sim <- array(NA, dim = dim_data)


  for (r in 1:n_rep) {
    print(r)
    X_train <- train_dist(n)
    eps_Y <- eps_dist(n)
    Y_train <- f_true(X_train) + eps_Y
    Y <- Y_train

    additive_pro <- progression::fit_additive(
      X = X_train,
      Y = Y_train,
      nb.iter = nb.iter,
      method = "llf_progression",
      pred_median = TRUE,
      plot = FALSE,
      # reg.fun = list(f_1 = reg_fun1, f_2 = reg_fun2),
      q = q,
      both.sided = TRUE,
      honesty = honesty_use,
      covariates = ncol(X_train)
    )

    rf <- regression_forest(
      X = X_train,
      Y = Y_train,
      num.trees = B0,
      honesty = honesty_use,
      min.node.size = 5
    )

    llf <- ll_regression_forest(
      X = X_train,
      Y = Y_train,
      enable.ll.split = TRUE,
      num.trees = B0,
      honesty = honesty_use,
      min.node.size = 5
    )

    en <- engression(X_train, Y_train, num_epochs = 1000)


    compute_error <- function(X_test) {
      Y_test <- f_true(X_test) + eps_dist(nrow(X_test))

      llf_pred <- predict(llf, X_test)$predictions
      forest_pred <- predict(rf, X_test)$predictions
      llf_prog_add_pred <- progression:::predict.custom_model(additive_pro, X_test)
      en_pred <- predict(en, X_test)

      pred_matrix <- cbind(
        forest_pred, llf_pred, en_pred, llf_prog_add_pred
      )

      RMSE <- apply(pred_matrix, MARGIN = 2, FUN = function(x) sqrt(mean((x - Y_test)^2)))
      MAE <- apply(pred_matrix, MARGIN = 2, FUN = function(x) mean(abs(x - Y_test)))
      MRAE <- apply(pred_matrix, MARGIN = 2, FUN = function(x) mean(pmin(abs(x - Y_test), abs(x / Y_test - 1))))

      return(cbind(RMSE, MAE, MRAE))
    }

    if (sim_setting == "polynomial") {
      X_test <- train_dist(n_test)
      res_sim[, 1, r, ] <- compute_error(X_test)

      X_test <- 3 * train_dist(n_test)
      res_sim[, 2, r, ] <- compute_error(X_test)

      X_test <- train_dist(n_test) + cbind(rep(0, times = n_test), rep(5, times = n_test))
      res_sim[, 3, r, ] <- compute_error(X_test)

      rho_test <- 0.8
      X_test <- 2 * qt(apply(rmvnorm(n, sigma = cbind(c(1, rho_test), c(rho_test, 1))), 2, unif), df = 4)
      res_sim[, 4, r, ] <- compute_error(X_test)
    }

    if (sim_setting == "Friedman") {
      X_test <- train_dist(n_test)
      res_sim[, 1, r, ] <- compute_error(X_test)

      X_test <- 3 * train_dist(n_test)
      res_sim[, 2, r, ] <- compute_error(X_test)

      X_test <- train_dist(n_test) + cbind(rep(1, times = n_test), rep(0, times = n_test), rep(0, times = n_test), rep(0, times = n_test), rep(0, times = n_test))
      res_sim[, 3, r, ] <- compute_error(X_test)


      X_test <- train_dist(n_test) + cbind(rep(0, times = n_test), rep(1, times = n_test), rep(0, times = n_test), rep(0, times = n_test), rep(0, times = n_test))
      res_sim[, 4, r, ] <- compute_error(X_test)
    }


  }
  save(res_sim, file=here("results", paste0("res_sim_run1_", sim_setting,".Rdata")))

}

#

load(here("results", paste0("res_sim_run1_polynomial.Rdata")))
box_data_poly = res_sim[,-1,,2]
load(here("results", paste0("res_sim_run1_Friedman.Rdata")))
box_data_Fried = res_sim[,-1,-9,2]



df1 <- as.data.frame(as.table(box_data_poly))
colnames(df1) <- c("Method", "Group", "Observation", "Value")
df1$Method <- factor(df1$Method, labels = c("RF","LLF","engression", "progression"))
df1$Group <- factor(df1$Group, labels = c("Var shift", "Mean shift", "Cov shift"))
df1 = df1 %>% mutate(h = "1")

df2 <- as.data.frame(as.table(box_data_Fried))
colnames(df2) <- c("Method", "Group", "Observation", "Value")
df2$Method <- factor(df2$Method, labels = c("RF","LLF","engression", "progression"))
df2$Group <- factor(df2$Group, labels = c("Var shift", "Mean shift X3", "Mean shift X4"))
df2 = df2 %>% mutate(h = "2")


df = bind_rows(df1,df2) %>% mutate(h = if_else(h=="1", "Polynomial", "Friedman"))
# df$Method = factor(df0$Method, levels = c("RF", "LLF", "engression", "progression"))
df$h = factor(df$h, levels=c("Polynomial", "Friedman"))

df1 = df1 %>% mutate(h = if_else(h=="1", "Polynomial", "Friedman"))


gg_box <- ggplot(df1, aes(x = Group, y = Value, fill = Method, color = Method)) +
  geom_boxplot(alpha = 0.8, size = 0.35, outlier.size = 0.7, width = 0.4) +
  facet_wrap(~h, , scales = "free_y") +
  xlab("") +
  ylab("RMSE") +
  scale_fill_manual(values = my_fill) +
  scale_colour_manual(values = my_col) +
  theme(
    legend.position = c(0.5, -0.23),
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    # strip.text.x = element_blank()
  )
gg_box

save_myplot(gg_box, plt_nm = here("num_results/plots", paste0("boxplot_MAE_poly.pdf")), width = 4, height = 4)


df2 = df2 %>% mutate(h = if_else(h=="1", "Polynomial", "Friedman"))


gg_box <- ggplot(df2, aes(x = Group, y = Value, fill = Method, color = Method)) +
  geom_boxplot(alpha = 0.8, size = 0.35, outlier.size = 0.7, width = 0.4) +
  facet_wrap(~h, , scales = "free_y") +
  xlab("") +
  ylab("RMSE") +
  scale_fill_manual(values = my_fill) +
  scale_colour_manual(values = my_col) +
  theme(
    legend.position = c(0.5, -0.23),
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    # strip.text.x = element_blank()
  )
gg_box

save_myplot(gg_box, plt_nm = here("num_results/plots", paste0("boxplot_MAE_fried.pdf")), width = 4, height = 4)


