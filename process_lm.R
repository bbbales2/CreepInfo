library(tidyverse)
library(rstan)
library(ggplot2)
library(readr)
library(rstanarm)
library(shinystan)
library(purrrlyr)

df = tbl_df(read_csv('/home/bbales2/CreepInfo/creep.csv')) %>%
  select(-iminf, -imaxf) %>%
  mutate(csv = paste(file, ".csv", sep = ""))

getSlope = function(row) {
  data = read_csv(row['csv'])
  
  fit = lm(stress ~ time, data = data)
  
  c(coefficients(fit)[['time']], sqrt(vcov(fit)[['time', 'time']]))
}

getSlope(as.matrix(df)[1,])

a = df %>% apply(1, getSlope) %>% t()
colnames(a)[1:2] <- c("slopes", "slope_std_err")

df = bind_cols(df, tbl_df(a))

df2 = df %>%
  select(thickness, stress, slopes, slope_std_err) %>%
  mutate(thickness_id = as.integer(as.factor(thickness))) %>%
  mutate(lmus = log(slopes))

df2 %>% print(n = 40)

data = list(L = nrow(df),
            T = max(df2$thickness_id),
            labels = df2$thickness_id,
            stress = df2$stress,
            thickness = df2$thickness,
            mus = df2$lmus)

fit = stan("/home/bbales2/CreepInfo/lumped.stan", data = data, cores = 4)

launch_shinystan(fit)

s = extract(fit)

colnames(s$muhat) <- 1:30

posterior = as.tibble(s$muhat[1:100,]) %>% gather("row", "lmus_calc", 1:30) %>% mutate(row = as.integer(row))

df3 = left_join(df2 %>% mutate(row = row_number()), posterior, by = "row")
df3 = df3 %>% mutate(thickness = factor(thickness))
df4 = df3 %>% group_by(thickness, stress) %>% summarize(lmus = max(lmus), lq = quantile(lmus_calc, 0.125), tq = quantile(lmus_calc, 0.875), lmus_calc = max(lmus_calc))

df3 %>% ggplot(aes(stress, exp(lmus_calc))) +
  geom_jitter(alpha = 0.2) +
  geom_errorbar(data = df4, aes(ymin = exp(lq), ymax = exp(tq)), color = "dodgerblue1", alpha = 0.8, size = 0.75) +
  geom_smooth(method = 'lm', formula = y ~ x, color = "dodgerblue1") +
  geom_point(data = df4, aes(stress, exp(lmus)), color = "red", shape = 17) +
  scale_y_log10() +
  facet_wrap(~ thickness, scales = "free", labeller = "label_both")

