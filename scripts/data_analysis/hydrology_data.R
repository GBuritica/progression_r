################################################
## Hydrology application
source('scripts/data_analysis/aux/applications.R') ## fct. to run all models
source('scripts/data_analysis/aux/stationarity_adjust.R') ## sts. adjust
source('scripts/plots/save_plot.R') ## saving plots
library(ggplot2)
library(egg)
library(latex2exp)
##
################################################
## load data
file        <- "data_analysis/hydrology_data/"
## Computing residuals (data after pre-treatment below)
## run to create data2.Rdata file
##precip_data      <- read.csv(paste0(file,"precip.csv")) ## Records per station
## discharges_data  <- read.csv(paste0(file,"Discharges_all.csv"))
## head(precip_data)
## stationarity adjust : season - trend
## d = 10, k = 20
#precip_adjusted     <- stationarity_adjust(precip_data,d = 10,k = 20)
#discharges_adjusted <- stationarity_adjust(discharges_data,d = 10,k = 20)
## extracts original data
#precip_records     <- precip_adjusted$data_original
#discharges_records <- discharges_adjusted$data_original
## extract adjusted data
#precip             <- precip_adjusted$data_adjusted
#discharges         <- discharges_adjusted$data_adjusted
#save( precip_adjusted, discharges_adjusted, precip_records, discharges_records, precip, discharges, file=paste0(data_dest_folder, '/hydrology_data/data2.Rdata') )
load(file = paste0('scripts/data_analysis/hydrology_data/data.Rdata') )
#######################################################
## Extracting non NA data - valid years
#err <- NULL

index_start <- which( discharges$Date == '1940-01-01')
index_end   <- which( discharges$Date == '2014-12-30') ## stop just before the August 22nd extreme event

##
ind_discharges <- index_start:index_end
##plots data
short_discharges_records <- discharges_records[ind_discharges, ]
g1 <- ggplot(data=short_discharges_records ) +
  geom_line(aes(x= Date, y=station_43)) + theme_bw() +
  labs(title = ' ',
       x = "Date",
       y = "Station 43")
g2 <- ggplot(data=short_discharges_records ) +
  geom_line(aes(x= Date, y=station_62)) + theme_bw() +
  labs(title = ' ',
       x = "Date",
       y = "Station 62")

#save_myplot(g1, paste0('scripts/hydrology_data/ts_records.pdf' ), width = 5, height =3)
#save_myplot(g2, paste0('scripts/hydrology_data/ts_records1.pdf' ), width = 5, height =3)

##plots residual data
short_discharges <- discharges[ind_discharges, ]
short_precip     <- precip[ind_precip, ]
g1 <- ggplot(short_discharges) +
  geom_line( aes(x= Date, y=station_43) ) + theme_bw() +
  labs(title = ' ',
       x = "Date",
       y = "Station 43")
g2 <- ggplot(short_discharges ) +
  geom_line(aes(x= Date, y=station_62)) + theme_bw() +
  labs(title = ' ',
       x = "Date",
       y = "Station 62")
#save_myplot(g1, paste0('scripts/hydrology_data/ts2.pdf' ), width = 5, height =3)
#save_myplot(g2, paste0('scripts/hydrology_data/ts22.pdf' ), width = 5, height =3)

n      <- length(short_discharges$Date)
data   <- data.frame(Date = short_discharges$Date[-c(1,2,3,4) ],
                     station_62 = short_discharges$station_62[-c(1,2,3,4) ],
                     #station_62_lag1  = short_discharges$station_62[ -n ], ## upstream station lag 1
                     station_43_lag3  = short_discharges$station_43[-c(1,n,n-1,n-2)], ## upstream station lag 1
                     station_43_lag2  = short_discharges$station_43[-c(1,2,n,n-1)], ## upstream station lag 1
                     station_43_lag1  = short_discharges$station_43[-c(1,2,3,n)] ## upstream station lag 1
)
p                   <- ncol(data)-2

######################################################
######################################################
## Regression model
## train index
date_start_train <- '1952-01-07' # date start train
date_train <- '1956-01-07' # date end train
## test index
date_start_test <- '2000-01-01'
date_end_test   <- '2014-01-01'


#err <- NULL # vector for saving the relative error test
{
set.seed(123)
## test index
index_start_test    <- which( data$Date == date_start_test)
index_end_test      <- which( data$Date == date_end_test)
index_start_train   <- which( data$Date == date_start_train)
index_end_train     <- which( data$Date == date_train)
ind_test            <- index_start_test:index_end_test
ind_train           <- 1:index_end_train

# Keep train vs. test clean without overlap
data$split <- NA_character_
data$split[ind_test]  <- "Test"
data$split[ind_train] <- "Train"

test_df  <- subset(data, split == "Test")
train_df <- subset(data, split == "Train")

g1 <- ggplot(mapping = aes(x = station_43_lag2, y = station_62)) +
  geom_point(
    data = test_df,
    aes(color = "Test"),
    alpha = 0.7,
    size = 1.8
  ) +
  geom_point(
    data = train_df,
    aes(color = "Train"),
    alpha = 0.7,
    size = 1.8
  ) +
  scale_color_manual(
    values = c("Test" = "gray60", "Train" = "black"),
    breaks = c("Train", "Test"),
    name = " "
  ) +
  xlim(quantile(train_df$station_43_lag2,.8),  150) +
  ylim(quantile(train_df$station_62,.8), 400) +
  labs(
    x = "Station 43 (lag 2)",
    y = "Station 62",
    title = " "
  ) +
  theme_bw()
g1

#save_myplot(g1, paste0('scripts/hydrology_data/train_test2.pdf' ), width = 5, height =3)
############################# Training
## training, test
X_train <- as.matrix(data[ind_train , 2 + (1:p)])
Y_train <- data$station_62[ind_train]
X_test  <- as.matrix(data[ind_test , 2 + (1:p)])
Y_test  <-data$station_62[ind_test]
## training all models
models     <- prediction_models(X_train,Y_train,X_test,Y_test, nb.iter= 1, cov = 3) #nb.iter= 50 for producing paper results
df_models0 <- data.frame(engression = models$engression$predictions,
                         progression = models$llf_additive$predictions,
                         forest = models$forest$predictions ,
                         llf = models$llf$predictions ,
                         true = Y_test,
                         X = X_test[,2])
df_models <- df_models0 # a copy for computing predictions at the real scale
head(df_models0)
##
index_start_test2    <- which( discharges_adjusted$data_season$Date == date_start_test)
index_end_test2      <- which( discharges_adjusted$data_season$Date == date_end_test)
ind_test2            <- index_start_test2:index_end_test2
ind_test2            <- ind_test2
##
index_start_test3    <- which( discharges_adjusted$data_trend$Date == date_start_test)
index_end_test3      <- which( discharges_adjusted$data_trend$Date == date_end_test)
ind_test3            <- index_start_test3:index_end_test3
ind_test3            <- ind_test3
### Estimates
df_models$engression    <- df_models0$engression + discharges_adjusted$data_season$station_62[ind_test2] + discharges_adjusted$data_trend$station_62[ind_test3]
df_models$progression   <- df_models0$progression + discharges_adjusted$data_season$station_62[ind_test2] + discharges_adjusted$data_trend$station_62[ind_test3]
df_models$forest        <- df_models0$forest + discharges_adjusted$data_season$station_62[ind_test2] + discharges_adjusted$data_trend$station_62[ind_test3]
df_models$llf           <- df_models0$llf + discharges_adjusted$data_season$station_62[ind_test2] + discharges_adjusted$data_trend$station_62[ind_test3]
df_models$true          <- df_models0$true + discharges_adjusted$data_season$station_62[ind_test2] +  discharges_adjusted$data_trend$station_62[ind_test3]
df_models$date          <- discharges_adjusted$data_season$Date[ind_test2]


## scatter plot
g1 <- ggplot(data=df_models) + geom_point(aes(x = true, y = engression), color = '#009E73', alpha = .6) + theme_bw() +
  geom_abline(intercept=0,slope=1)+
  labs(title = 'Engression',
       y =TeX("Predictions ($m^3/s$)"), x = TeX("True values ($m^3/s$)")) + xlim(0,800)  + ylim(0,800)
g2 <- ggplot(data=df_models) + geom_point(aes(x = true, y = progression), color = '#E69F00', alpha = .6) + theme_bw() +
  geom_abline(intercept=0,slope=1)+
  labs(title = 'Progression',
       y =TeX("Predictions ($m^3/s$)"), x = TeX("True values ($m^3/s$)")) + xlim(0,800)  + ylim(0,800)
g3 <- ggplot(data=df_models) + geom_point(aes(x = true, y = forest), color = '#D55E00', alpha = .6) + theme_bw() +
  geom_abline(intercept=0,slope=1)+
  labs(title = 'Random Forest',
       y =TeX("Predictions ($m^3/s$)"), x = TeX("True values ($m^3/s$)")) + xlim(0,800)  + ylim(0,800)
g4 <- ggplot(data=df_models) + geom_point(aes(x = true, y = llf), color = '#0072B2', alpha = .6) + theme_bw() +
  geom_abline(intercept=0,slope=1)+
  labs(title = 'LLF',
       y =TeX("Predictions ($m^3/s$)"), x = TeX("True values ($m^3/s$)")) + xlim(0,800)  + ylim(0,800)
ggarrange(g1,g2,g3,g4,ncol=2)

#save_myplot(g1, paste0('scripts/hydrology_data/scatter_engression_all2', date_train, '.pdf' ), width = 5, height =4)
#save_myplot(g2, paste0('scripts/hydrology_data/scatter_progression_all2_rf2', date_train, '.pdf' ), width = 5, height =4)
#save_myplot(g3, paste0('scripts/hydrology_data/scatter_forest_all2', date_train, '.pdf' ), width = 5, height =4)
#save_myplot(g4, paste0('scripts/hydrology_data/scatter_llf_all2', date_train, '.pdf' ), width = 5, height =4)

## plot predictions
ind <- which(df_models$date=='2002-01-01'):which(df_models$date=='2008-01-01')
df_ind <- df_models[ind,]
df_ind$date <- data$Date[ind_test[ind] ]
g1 <- ggplot(df_ind, aes(x = date, y = engression)  ) +
  geom_point( aes(x =date, y = true) , color = 'black', alpha = .5, size=.2) +
  geom_line(color= '#009E73',alpha = .7) + ylim(30,700)+
  ylab(expression(paste("Discharge (", m^3, "/s)"))) + xlab('Date')+
  labs(title = "Engression")+
  theme_bw()
g3 <- ggplot(df_ind, aes(x = date, y = llf)  ) +
  geom_point( aes(x =date, y = true) , color = 'black', alpha = .5, size=.2) +
  geom_line(color= '#0072B2',alpha = .7) + ylim(30,700)+
  ylab(expression(paste("Discharge (", m^3, "/s)"))) + xlab('Date')+
  labs(title = "LLF")+
  theme_bw()
g4 <- ggplot(df_ind, aes(x = date, y = forest)  ) +
  geom_point( aes(x =date, y = true) , color = 'black', alpha = .5, size=.2) +
  geom_line(color= '#D55E00',alpha = .7) + ylim(30,700)+
  ylab(expression(paste("Discharge (", m^3, "/s)"))) + xlab('Date')+
  labs(title = "Random Forest")+
  theme_bw()
g2 <- ggplot(df_ind, aes(x = date, y = progression)  ) +
  geom_point( aes(x =date, y = true) , color = 'black', alpha = .5, size=.2) +
  geom_line(color= '#E69F00',alpha = .7) + ylim(30,700)+
  ylab(expression(paste("Discharge (", m^3, "/s)"))) + xlab('Date')+
  labs(title = "Progression")+
  theme_bw()

grid.arrange(g1,g2,g3,g4,nrow=2)

#save_myplot(g1, paste0(figure_dest_folder, '/hydrology_data/ts_engression_all2title', date_train, '.pdf' ), width = 4, height =3)
#save_myplot(g2, paste0(figure_dest_folder, '/hydrology_data/ts_progression_all2title', date_train, '.pdf' ), width = 4, height =3)
#save_myplot(g4, paste0(figure_dest_folder, '/hydrology_data/ts_forest_all2title', date_train, '.pdf' ), width = 4, height =3)
#save_myplot(g3, paste0(figure_dest_folder, '/hydrology_data/ts_llf_all2title', date_train, '.pdf' ), width = 4, height =3)

## For each test period we compute the relative error
## relative error
true <- df_models$true
sor  <- sort(true)
for(i in 1:length(sor)){
  ind <- which(true == sor[i] )
  err <- rbind(err, c(apply(df_models[ind,1:5],2, function(x) mean(abs(x- df_models$true[ind])/df_models$true[ind] ) ), sor[i] ))
}

}



#save(err, true, file = paste0('scripts/hydrology_data/err', date_train, '.Rdata' ))
#load(file = paste0('scripts/data_analysis/hydrology_data/err', '1960-01-07', '.Rdata' ))
err <- as.data.frame(err)

library(dplyr)

median_traj <- err %>%
  group_by(V6) %>%
  summarise(
    engressionqu = quantile(engression, probs= .9, na.rm = TRUE),
    progressionqu = quantile(progression, probs= .9,na.rm = TRUE),
    llfqu = quantile(llf,  probs=.9, na.rm = TRUE),
    forestqu = quantile(forest, probs= .9, na.rm = TRUE),
    engressionql = quantile(engression,  probs=.1, na.rm = TRUE),
    progressionql = quantile(progression,  probs=.1,na.rm = TRUE),
    llfql = quantile(llf, probs= .1, na.rm = TRUE),
    forestql = quantile(forest,  probs= .1, na.rm = TRUE),
    engression = mean(engression,na.rm = TRUE),
    progression = mean(progression,na.rm = TRUE),
    llf = mean(llf, na.rm = TRUE),
    forest = mean(forest, na.rm = TRUE)
  ) %>%
  arrange(V6) %>%
  ungroup()


g1 <- ggplot(median_traj, aes(x=V6)) +
  geom_point(aes(y=engression), color= '#009E73' ,alpha=.7, cex=.7) +
  xlim(250,560)+
  ylim(0,1)+
  geom_point(aes(y=forest), color = '#D55E00',alpha=.7, cex=.7) +
  geom_point(aes(y=llf), color = '#0072B2',alpha=.7, cex=.7) +
  geom_point(aes(y=progression), color = '#E69F00',alpha=.7, cex=.7) +
  geom_line(aes(y=engression), color= '#009E73' ,alpha=.2) +
  geom_line(aes(y=forest), color = '#D55E00',alpha=.2) +
  geom_line(aes(y=progression), color = '#E69F00',alpha=.2) +
  geom_line(aes(y=llf), color = '#0072B2',  alpha=.2) +
  geom_segment(aes(x = V6[5045]-1.5, y =  engressionqu[5045], xend = V6[5045]-1.5, yend =  engressionql[5045]),color='#009E73',  alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5045], y =  forestqu[5045], xend = V6[5045], yend =  forestql[5045]),color='#D55E00',  alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5045]+.5, y =  llfqu[5045], xend = V6[5045]+.5, yend =  llfql[5045]),color='#0072B2',alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5045]+2.5, y =  progressionqu[5045], xend = V6[5045]+2.5, yend =  progressionql[5045]),color='#E69F00',  alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5040]-1.5, y =  engressionqu[5040], xend = V6[5040]-1.5, yend =  engressionql[5040]),color='#009E73',  alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5040], y =  forestqu[5040], xend = V6[5040], yend =  forestql[5040]),color='#D55E00',  alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5040]+.5, y =  llfqu[5040], xend = V6[5040]+.5, yend =  llfql[5040]),color='#0072B2',alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5040]+2.5, y =  progressionqu[5040], xend = V6[5040]+2.5, yend =  progressionql[5040]),color='#E69F00',  alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5033]-1, y =  engressionqu[5033], xend = V6[5033]-1, yend =  engressionql[5030]),color='#009E73', alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5033], y =  forestqu[5033], xend = V6[5033], yend =  forestql[5030]),color='#D55E00',  alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5033]+.5, y =  llfqu[5033], xend = V6[5033]+.5, yend =  llfql[5030]),color='#0072B2',alpha=.3, linewidth = .5)+
  geom_segment(aes(x = V6[5033]+2.3, y =  progressionqu[5033], xend = V6[5033]+2.3, yend =  progressionql[5030]),color='#E69F00',  alpha=.3, linewidth = .5)+
  ylab('Relative test error') +
  xlab(TeX("Threshold $u$ ($m^3/s$)"))+
  theme_bw()

g1

#save_myplot(g1, paste0('scripts/data_analysis/ts_llf_all2title', date_train, '.pdf' ), width = 4, height =3)

