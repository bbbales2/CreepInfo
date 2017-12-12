library(tidyverse)
library(rstan)
library(ggplot2)
library(readr)
library(rstanarm)
library(shinystan)
library(purrrlyr)
library(plotly)
library(ggthemes)

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

#launch_shinystan(fit)

s = extract(fit)

dfp = left_join(df2 %>% mutate(idx = 1),
                as_tibble(list(aa = s$a, bb = s$b, cc = s$c)) %>% mutate(idx = 1)) %>%
  select(-idx) %>%
  group_by(thickness) %>%
  sample_n(1000) %>%
  ungroup()

dfp %>% ggplot(aes(exp(cc) * (stress^aa), slopes / ((1 / thickness)^bb))) +
  geom_point(aes(color = factor(thickness)), alpha = 0.25) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() + guides(colour = guide_legend(override.aes = list(size = 5)))

colnames(s$muhat) <- 1:30
colnames(s$mumu) <- 1:30

posterior = bind_cols(as.tibble(s$muhat[,]) %>%
                        gather("row", "lmus_calc", 1:30) %>%
                        mutate(row = as.integer(row)),
                      as.tibble(s$mumu[,]) %>%
                        gather("row", "mumus_calc", 1:30) %>%
                        select(-row))

df3 = left_join(df2 %>% mutate(row = row_number()), posterior, by = "row")
df3 = df3 %>% mutate(thickness = factor(thickness))
df4 = df3 %>% group_by(thickness, stress) %>%
  summarize(mean_mumus_calc = mean(mumus_calc),
            mean_mumus_m2sd = mean_mumus_calc - 2.0 * sd(mumus_calc),
            mean_mumus_p2sd = mean_mumus_calc + 2.0 * sd(mumus_calc),
            mean_lmus_calc = mean(lmus_calc),
            mean_lmus_m2sd = mean_lmus_calc - 2.0 * sd(lmus_calc),
            mean_lmus_p2sd = mean_lmus_calc + 2.0 * sd(lmus_calc),
            lmus_calc = NA) %>%
  right_join(df2 %>% mutate(thickness = factor(thickness)), by = c("thickness", "stress"))

df3 %>% group_by(thickness, stress) %>%
  summarize(sd_abc = sd(mumus_calc) / log(10), sd_abc_sigma = sd(lmus_calc) / log(10)) %>%
  print(n = 40) %>%
  group_by(thickness)

df3 %>% group_by(thickness, stress) %>%
  summarize(sd_abc = sd(mumus_calc) / log(10), sd_abc_sigma = sd(lmus_calc) / log(10)) %>%
  group_by(thickness) %>%
  summarize(mu_sd_abc = mean(sd_abc), mean_sd_abc_sigma = mean(sd_abc_sigma))
  
df3 %>% group_by(thickness, stress) %>% sample_n(100) %>% ggplot(aes(stress, exp(lmus_calc))) +
  geom_jitter(alpha = 0.2) +
  #geom_errorbar(data = df4, aes(ymin = exp(m2sd), ymax = exp(p2sd)), color = "dodgerblue1", alpha = 0.8, size = 0.75) +
  geom_line(data = df4, aes(stress, exp(mean_mumus_calc)), color = "dodgerblue1") +
  geom_ribbon(data = df4, aes(stress, ymin = exp(mean_lmus_m2sd), ymax = exp(mean_lmus_p2sd)), fill = "grey", alpha = 0.4) +
  geom_ribbon(data = df4, aes(stress, ymin = exp(mean_mumus_m2sd), ymax = exp(mean_mumus_p2sd)), fill = "dodgerblue1", alpha = 0.2) +
  geom_point(data = df4, aes(stress, exp(lmus)), color = "violetred2") +
  scale_y_log10() +
  facet_wrap(~ thickness, labeller = "label_both") + theme_excel()# + scale_color_excel()

df4 %>% mutate(ln_slopes = lmus) %>% select(-lmus_calc, -thickness_id, -lmus) %>% print(n = 40)

colnames(s$sigma) <- 1:4

df4 %>% group_by(thickness, stress) %>%

as.tibble(s$sigma[,]) %>%
  gather("row", "sigmas", 1:4) %>%
  mutate(row = as.integer(row)) %>%
  left_join()

