library(tidyverse)
library(rstan)
library(ggplot2)
library(readr)
library(shinystan)
#library(purrrlyr)

# Read in all the list of files to process
df = tbl_df(read_csv('creep2.csv'))

getSlope = function(filename) {
  data = read_csv(filename)
  
  fit = lm(stress ~ time, data = data)
  
  c(coefficients(fit)[['time']], sqrt(vcov(fit)[['time', 'time']]))
}

# For each csv file, read in the data, get a linear fit, and put it in a dataframe (df2)
a = df %>%
  pull(csv) %>%
  map(getSlope) %>%
  unlist() %>%
  matrix(ncol = 2, byrow = TRUE)
colnames(a)[1:2] <- c("slopes", "slope_std_err")

df2 = bind_cols(df, tbl_df(a)) %>%
  select(thickness, stress, slopes, slope_std_err) %>%
  mutate(thickness_id = as.integer(as.factor(thickness))) %>%
  mutate(lmus = log(slopes))

df2 = df2 %>% left_join(df2 %>% group_by(thickness) %>%
                          summarize(mean_log_stress = mean(log(stress))))

# This prints the slopes and errors in slopes for all combos of thicknesses and stress
df2 %>% print(n = 40)

data = list(L = nrow(df2),
            T = max(df2$thickness_id),
            labels = df2$thickness_id,
            log_stress = log(df2$stress),# - mean(log(df2$stress)),#df2$mean_log_stress,
            log_inv_thickness = log(1 / df2$thickness),# - mean(log(1 / c(65, 200, 500, 2000))),
            mus = df2$lmus,
            log_inv_thicknesses = log(1 / c(65, 200, 500, 2000)))# - mean(log(1 / c(65, 200, 500, 2000))))

# Run the LMP fit
fit = stan("/home/bbales2/CreepInfo/lumped.stan", data = data, cores = 4)
fit2 = stan("/home/bbales2/CreepInfo/hierarchical.stan", data = data, cores = 4)
fit3 = stan("/home/bbales2/CreepInfo/full_hierarchical.stan", data = data, cores = 4,
            control = list(max_treedepth = 12))

launch_shinystan(fit)
launch_shinystan(fit2)
launch_shinystan(fit3)

s = rstan::extract(fit)

gquants = function(a) {
  a %>%
    as.tibble %>%
    setNames(1:nrow(df2)) %>%
    gather(row, y) %>%
    mutate(row = as.integer(row)) %>%
    left_join(df2 %>% mutate(row = row_number()), by = "row") %>%
    group_by(thickness, stress) %>%
    summarize(y5 = quantile(y, 0.05),
              y95 = quantile(y, 0.95))
}

s$sdo_hat %>%
  gquants %>%
  ggplot(aes(stress)) +
  geom_ribbon(aes(ymin = y5, ymax = y95, group = thickness, fill = as.factor(thickness)), alpha = 0.25) +
  geom_errorbar(data = s$sdo %>% gquants, aes(ymin = y5, ymax = y95, group = thickness, color = as.factor(thickness))) +
  geom_point(data = df2, aes(stress, lmus - mean(s$p) * log(1 / df2$thickness)))

s$lso_hat %>%
  gquants %>%
  mutate(measured_slopes = TRUE) %>%
  bind_rows(s$lso %>% gquants %>% mutate(measured_slopes = FALSE)) %>%
  ggplot(aes(log(1 / thickness))) +
  geom_errorbar(aes(ymin = y5, ymax = y95, group = stress, color = as.factor(thickness)), alpha = 0.25, size = 1) +
  geom_point(data = df2, aes(log(1 / thickness), lmus - mean(s$n) * log(df2$stress))) +
  facet_grid(. ~ measured_slopes, labeller = "label_both")

# Get all the data out and plot it!
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

colnames(s$muhat) <- 1:nrow(df2)
colnames(s$mumu) <- 1:nrow(df2)

posterior = bind_cols(as.tibble(s$muhat[,]) %>%
                        gather("row", "lmus_calc", 1:nrow(df2)) %>%
                        mutate(row = as.integer(row)),
                      as.tibble(s$mumu[,]) %>%
                        gather("row", "mumus_calc", 1:nrow(df2)) %>%
                        select(-row))

df3 = left_join(df2 %>% mutate(row = row_number()), posterior, by = "row")
df3 = df3 %>% mutate(thickness = factor(thickness))
df4 = df3 %>% group_by(thickness, stress) %>%
  summarize(mean_mumus_calc = mean(mumus_calc),
            mean_mumus_025 = quantile(mumus_calc, 0.025),#mean_mumus_calc - 2.0 * sd(mumus_calc),
            mean_mumus_975 = quantile(mumus_calc, 0.975),
            mean_lmus_calc = mean(lmus_calc),
            mean_lmus_025 = quantile(lmus_calc, 0.025),
            mean_lmus_975 = quantile(lmus_calc, 0.975),
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
  #geom_jitter(alpha = 0.2) +
  #geom_errorbar(data = df4, aes(ymin = exp(025), ymax = exp(975)), color = "dodgerblue1", alpha = 0.8, size = 0.75) +
  geom_line(data = df4, aes(stress, exp(mean_mumus_calc)), color = "dodgerblue1") +
  geom_ribbon(data = df4, aes(stress, ymin = exp(mean_lmus_025), ymax = exp(mean_lmus_975)), fill = "grey", alpha = 0.4) +
  geom_ribbon(data = df4, aes(stress, ymin = exp(mean_mumus_025), ymax = exp(mean_mumus_975)), fill = "dodgerblue1", alpha = 0.2) +
  geom_point(data = df4, aes(stress, exp(lmus)), color = "violetred2") +
  scale_y_log10() +
  facet_wrap(~ thickness, labeller = "label_both")# + scale_color_excel()

df4 %>% mutate(ln_slopes = lmus) %>% select(-lmus_calc, -thickness_id, -lmus) %>% print(n = 40)

udf = s$uncertainty %>%
  as.tibble %>%
  setNames(1:ncol(.)) %>%
  gather(row, uncertainty) %>%
  mutate(row = as.integer(row),
         uncertainty = uncertainty / log(10)) %>%
  left_join(df2 %>% mutate(row = row_number()), posterior, by = "row")

udf %>%
  ggplot(aes(stress, uncertainty)) +
  geom_jitter(alpha = 0.01) +
  facet_wrap(~ thickness)

(udf %>% group_by(thickness, stress) %>%
    summarize(m = mean(uncertainty),
              sd = sd(uncertainty),
              qm1 = quantile(uncertainty, probs = c(0.159)),
              qp1 = quantile(uncertainty, probs = c(0.841))) %>% print(n = nrow(.))) %>%
  ggplot(aes(stress, m)) +
  geom_errorbar(aes(stress, ymin = m - sd, ymax = m + sd), colour = "red") +
  geom_errorbar(aes(stress, ymin = qm1, ymax = qp1), colour = "green") +
  geom_point() +
  facet_wrap(~ thickness, labeller = "label_both") +
  ggtitle("mean of log(pred) - log(actual) black dot\ngreen estimated +- 1 quantile\nred estimated +- 1 quantile assuming difference is normally distributed")
