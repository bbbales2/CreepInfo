library(tidyverse)
library(rstan)
library(ggplot2)
library(readr)
library(rstanarm)
library(shinystan)
library(purrlyr)

df = tbl_df(read_csv('/home/bbales2/CreepInfo/creep.csv')) %>%
  select(-iminf, -imaxf) %>%
  mutate(csv = paste(file, ".csv", sep = ""))

df %>% apply(1, test)

df %>% mutate(asdf = read_csv(csv))

getSlope = function(row) {
  data = read_csv(row['csv'])
  
  data$time = data$time
  
  fit = stan_glm(stress ~ time, data = data, cores = 4)
  
  draws = as.data.frame(fit)
  
  colnames(draws)[1] <- c("intercept")
  
  plot = data %>% ggplot(aes(x = time, y = stress)) + 
    geom_abline(data = draws, aes(intercept = intercept, slope = time), 
                color = "skyblue", size = 0.2, alpha = 0.1) +
    geom_point() +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], 
                color = "red", size = 1) +
    ggtitle(row['csv'])
  
  list(mean(draws$time), sd(draws$time), mean(draws$sigma), fit, plot)
}

getSlope(as.matrix(df)[29,])

a = df %>% apply(1, getSlope) %>% t()
colnames(a)[1:3] <- c("slopes", "slope_std_err", "sigma")

df = bind_cols(df, tbl_df(a))

df %>% select(thickness, stress, slopes, slope_std_err, sigma) %>% print(n = 40)
