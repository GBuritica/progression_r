## script for producing Figure 1 from the paper
## progression : regression extrapolation principle
figure_dest_folder <- 'results/plots'
data_dest_folder <- 'results/data'
library(ggplot2)
library(ggpubr)
library(progression)
library(grf)
library(patchwork)

source('scripts/plots/save_plot.R') ##

set.seed(145)


n=300
sd_eps = 6
df_X = 6
f = function(x) .5*x^3 + x^2 - x
X = rt(n, df = df_X) #2*(rpareto(n, shape=df_X) -rpareto(n, shape=df_X)) #
eps =  rnorm(n, sd = sd_eps)
Y = f(X) + eps


#### population level approximations
#library(VGAM)  # for laplace distribution
#library(skewsamp) # for inverse of ecdf
#n_true=50000
#X_true = rt(n_true, df = df_X)
#eps_true =  rnorm(n_true, sd = sd_eps)
#Y_true = f(X_true) + eps_true
#ecdf_Y = ecdf(Y_true) #approximation of true cdf


#X_LP = qlaplace(pt(X, df=df_X))
#Y_LP = qlaplace(ecdf_Y(Y))
# HT_function = function(x) x

## GPD fit
q <- 0.08
B0 <- 2000

trans_X <- prepare_laplace_transformations(data = X, p = 1-q)
trans_Y <- prepare_laplace_transformations(data = Y, p = 1-q)
X_LP <- trans_X$data_to_laplace(X)
Y_LP <- trans_Y$data_to_laplace(Y)


uX <- quantile(X, p = 1-q)
uX_LP <- quantile(X_LP, p = 1-q)

## Fit at locally linear random forest
llf_LP <- ll_regression_forest(
  X = as.matrix(X_LP),
  Y = Y_LP,
  enable.ll.split = TRUE,
  num.trees = B0,
  honesty = FALSE,
  min.node.size = 10)


HT_fct = function(x) predict(llf_LP, as.matrix(x))$predictions

library(dplyr)
dat_orig <- tibble(X = X, Y=Y)

dat_LP <- tibble(X = X_LP, Y=Y_LP)
dat_LP <- dat_LP  %>%  mutate(HT = if_else(X > uX_LP, HT_fct(X), NaN))


dat_orig <- dat_orig  %>%  mutate(HT = if_else(X > uX, trans_Y$laplace_to_data(HT_fct(X_LP)) , NaN), f_true = f(X))
HT <- dat_orig$HT
f_true <- dat_orig$f_true

gg2 <-  ggplot(dat_orig) +
  geom_point(aes(x = X, y = Y), alpha = 0.2) +
  scale_color_manual(values = c( 'orange',"#0072B2"))+
  geom_line(aes(x = X, y = f_true),  color='orange') +
  coord_cartesian(xlim = c(min(X),max(X)),
                  ylim=c(min(Y, HT[!is.na(HT)], f_true),
                         max(Y, HT[!is.na(HT)], f_true))) +
  # theme_classic()+
  xlab(latex2exp::TeX('Predictor $X$' )) +
  ylab(latex2exp::TeX("Response $Y$")) +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.1),  # Cuadrícula principal
    panel.grid.minor = element_line(color = "gray", size = 0.1)  # Cuadrícula secundaria
  )
gg2

save_myplot(
  plt = gg2,
  plt_nm = paste0(figure_dest_folder, "plot_HT_1", ".pdf"),
  width = 3,
  height = 3
)

gg2a <-  ggplot(dat_orig) +
  geom_point(aes(x = X, y = Y), alpha = 0.2) +
  scale_color_manual(values = c( 'orange',"#0072B2"))+
  geom_line(aes(x = X, y = f_true),  color='orange') +
  annotate(geom="text", x=uX, y=-30, label=latex2exp::TeX('$u_X$'))  +
  geom_vline(xintercept=uX, linetype="dashed", color = 'darkgray') +
  # theme_classic()+
  xlab(latex2exp::TeX('Predictor $X$' )) +
  ylab(latex2exp::TeX("Response $Y$")) +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.1),  # Cuadrícula principal
    panel.grid.minor = element_line(color = "gray", size = 0.1)  # Cuadrícula secundaria
  )
gg2a

save_myplot(
  plt = gg2a,
  plt_nm = paste0(figure_dest_folder, "plot_HT_1a", ".pdf"),
  width = 3,
  height = 3
)


gg3 <-  ggplot(dat_LP) +
  geom_point(aes(x = X, y = Y), alpha = 0.2) +
  geom_line(aes(x = X, y = HT),  color="#0072B2") +
  geom_vline(xintercept=uX_LP, linetype="dashed", color='darkgray')+
  annotate(geom="text", x=uX_LP, y=-10.8, label=latex2exp::TeX('$u_X^*$'))  +
  coord_cartesian(xlim = c(min(X_LP),max(X_LP)),
                  ylim=c(min(Y_LP),max(Y_LP)), clip = "off") +
  # theme_classic()+
  xlab(latex2exp::TeX('Predictor $X^*$' )) +
  ylab(latex2exp::TeX("Response $Y^*$")) +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.1),  # Cuadrícula principal
    panel.grid.minor = element_line(color = "gray", size = 0.1)  # Cuadrícula secundaria
  )
gg3


save_myplot(
  plt = gg3,
  plt_nm = paste0(figure_dest_folder, "plot_HT_2", ".pdf"),
  width = 3,
  height = 3
)


gg4 <-  ggplot(dat_orig) +
  geom_point(aes(x = X, y = Y), alpha = 0.2) +
  scale_color_manual(values = c( 'orange',"#0072B2"))+
  geom_line(aes(x = X, y = f_true),  color='orange') +
  geom_line(aes(x = X, y = HT),  color="#0072B2") +
  annotate(geom="text", x=uX, y=-40, label=latex2exp::TeX('$u_X$'))  +      geom_vline(xintercept=uX, linetype="dashed", color = 'darkgray') +
  coord_cartesian(xlim = c(min(X),max(X)),
                  ylim=c(min(Y, HT[!is.na(HT)], f_true),max(Y, HT[!is.na(HT)], f_true)), clip = "off") +
  # theme_classic()+
  xlab(latex2exp::TeX('Predictor $X$' )) +
  ylab(latex2exp::TeX("Response $Y$")) +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.1),  # Cuadrícula principal
    panel.grid.minor = element_line(color = "gray", size = 0.1)  # Cuadrícula secundaria
  )
gg4

save_myplot(
  plt = gg4,
  plt_nm = paste0(figure_dest_folder, "plot_HT_3", ".pdf"),
  width = 3,
  height = 3
)



######################################
### Population level
### light tailed plot

set.seed(1234)
n=300
sd_eps = 2
scaleX = 7
f = function(x)  progression:::spline5(x)#x+0.1*x^2#
alpha = 3#2#
X = rweibull(n, shape= 2 * alpha , scale = scaleX * alpha) - rweibull(n, shape= 2 * alpha , scale = scaleX * alpha)
eps =  rnorm(n, sd = sd_eps)
Y = f(X) + eps
#curve(f,from = -15,to=15)
#points(X,Y)

# n_true=50000
# X_true = rweibull(n_true, shape= 2 * alpha , scale = scaleX * alpha) - rweibull(n_true, shape= 2 * alpha , scale = scaleX * alpha)
# eps_true =  rnorm(n_true, sd = sd_eps)
# Y_true = f(X_true) + eps_true
# ecdf_Y = ecdf(Y_true) #approximation of true cdf
# ecdf_X = ecdf(X_true)
#
# X_LP = qlaplace(ecdf_X(X))
# Y_LP = qlaplace(ecdf_Y(Y))

# par <- lm(Y_LP[X_LP>uX_LP]~X_LP[X_LP>uX_LP])
# par
# HT_fct = function(x) x * par$coefficients[2] + par$coefficients[1]


q <- 0.08
B0 <- 2000


trans_X <- prepare_laplace_transformations(data = X, p = 1-q)
trans_Y <- prepare_laplace_transformations(data = Y, p = 1-q)
X_LP <- trans_X$data_to_laplace(X)
Y_LP <- trans_Y$data_to_laplace(Y)
HT_fct = function(x) predict(llf_LP, as.matrix(x))$predictions

uX <- quantile(X, p = 1-p)
uX_LP <- quantile(X_LP, p = 1-p)

## Fit at locally linear random forest
llf_LP <- ll_regression_forest(
  X = as.matrix(X_LP),
  Y = Y_LP,
  enable.ll.split = TRUE,
  num.trees = B0,
  honesty = F,
  min.node.size = 10)


dat_orig <- tibble(X = X, Y=Y)


dat_LP <- tibble(X = X_LP, Y=Y_LP)
dat_LP <- dat_LP  %>%  mutate(HT = if_else(X > uX_LP, HT_fct(X), NaN))


dat_orig <- dat_orig  %>%  mutate(HT = if_else(X > uX, trans_Y$laplace_to_data(HT_fct(X_LP)) , NaN), f_true = f(X))
HT <- dat_orig$HT
f_true <- dat_orig$f_true



gg2 <-  ggplot(dat_orig) +
  geom_point(aes(x = X, y = Y), alpha = 0.2) +
  scale_color_manual(values = c( 'orange',"#0072B2"))+
  geom_line(aes(x = X, y = f_true),  color='orange') +
  coord_cartesian(xlim = c(min(X),max(X)),
                  ylim=c(min(Y, HT[!is.na(HT)], f_true),max(Y, HT[!is.na(HT)], f_true)), clip="off") +
  xlab(latex2exp::TeX('Predictor $X$' )) +
  ylab(latex2exp::TeX("Response $Y$")) +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.1),  # Cuadrícula principal
    panel.grid.minor = element_line(color = "gray", size = 0.1)  # Cuadrícula secundaria
  )
gg2

save_myplot(
  plt = gg2,
  plt_nm = paste0(figure_dest_folder, "1plot_HT_1", ".pdf"),
  width = 3,
  height = 3
)

gg2a <-  ggplot(dat_orig) +
  geom_point(aes(x = X, y = Y), alpha = 0.2) +
  scale_color_manual(values = c( 'orange',"#0072B2"))+
  geom_line(aes(x = X, y = f_true),  color='orange') +
  annotate(geom="text", x=uX, y=-25, label=latex2exp::TeX('$u_X$'))  +      geom_vline(xintercept=uX, linetype="dashed", color = 'darkgray') +
  coord_cartesian(xlim = c(min(X),max(X)),
                  ylim=c(min(Y, HT[!is.na(HT)], f_true),max(Y, HT[!is.na(HT)], f_true)), clip="off") +
  xlab(latex2exp::TeX('Predictor $X$' )) +
  ylab(latex2exp::TeX("Response $Y$")) +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.1),  # Cuadrícula principal
    panel.grid.minor = element_line(color = "gray", size = 0.1)  # Cuadrícula secundaria
  )
gg2a

save_myplot(
  plt = gg2a,
  plt_nm = paste0(figure_dest_folder, "1plot_HT_1a", ".pdf"),
  width = 3,
  height = 3
)


gg3 <-  ggplot(dat_LP) +
  geom_point(aes(x = X, y = Y), alpha = 0.2) +
  geom_line(aes(x = X, y = HT),  color="#0072B2") +
  geom_vline(xintercept=uX_LP, linetype="dashed", color='darkgray')+
  annotate(geom="text", x=uX_LP, y=-6.4, label=latex2exp::TeX('$u_X^*$'))  +
  coord_cartesian(xlim = c(min(X_LP),max(X_LP)),
                  ylim=c(min(Y_LP),max(Y_LP)), clip="off") +
  xlab(latex2exp::TeX('Predictor $X$' )) +
  ylab(latex2exp::TeX("Response $Y$")) +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.1),  # Cuadrícula principal
    panel.grid.minor = element_line(color = "gray", size = 0.1)  # Cuadrícula secundaria
  )
gg3


save_myplot(
  plt = gg3,
  plt_nm = paste0(figure_dest_folder, "1plot_HT_2", ".pdf"),
  width = 3,
  height = 3
)


gg4 <-  ggplot(dat_orig) +
  geom_point(aes(x = X, y = Y), alpha = 0.2) +
  scale_color_manual(values = c( 'orange',"#0072B2"))+
  geom_line(aes(x = X, y = f_true),  color='orange') +
  geom_line(aes(x = X, y = HT),  color="#0072B2") +
  annotate(geom="text", x=uX, y=-25, label=latex2exp::TeX('$\\tau$'))  +      geom_vline(xintercept=uX, linetype="dashed", color = 'darkgray') +
  coord_cartesian(xlim = c(min(X),max(X)),
                  ylim=c(min(Y, HT[!is.na(HT)], f_true),max(Y, HT[!is.na(HT)], f_true)), clip="off") +
  theme_classic()+
  xlab(latex2exp::TeX('Predictor $X$' )) +
  ylab(latex2exp::TeX("Response $Y$")) +
  theme(
    panel.grid.major = element_line(color = "gray", size = 0.1),  # Cuadrícula principal
    panel.grid.minor = element_line(color = "gray", size = 0.1)  # Cuadrícula secundaria
  )
gg4

save_myplot(
  plt = gg4,
  plt_nm = paste0(figure_dest_folder, "1plot_HT_3", ".pdf"),
  width = 3,
  height = 3
)


########################################################
########################################################


