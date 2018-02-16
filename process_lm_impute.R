library(tidyverse)
library(rstan)
library(ggplot2)
library(readr)

# Read in all the list of files to process
df = tbl_df(read_csv('creep2.csv'))

getSlope = function(filename) {
  data = read_csv(filename)
  
  fit = lm(stress ~ time, data = data)
  
  c(coefficients(fit)[['time']], sqrt(vcov(fit)[['time', 'time']]))
}

getSlopes = function(N, filename) {
  o = getSlope(filename)
  
  rnorm(N, o[1], o[2])
}

N = 100
# For each csv file, read in the data, get a linear fit, and put it in a dataframe (df2)
a = df %>%
  pull(csv) %>%
  map(~ getSlopes(N, .)) %>%
  unlist %>%
  matrix(ncol = N, byrow = TRUE)

fits = list()
for(n in 1:N) {
  df2 = bind_cols(df, tbl_df(a[, n]) %>% setNames("slopes")) %>%
    select(thickness, stress, slopes) %>%
    mutate(thickness_id = as.integer(as.factor(thickness))) %>%
    mutate(lmus = log(slopes))

  data = list(L = nrow(df2),
              T = max(df2$thickness_id),
              labels = df2$thickness_id,
              stress = df2$stress,
              thickness = df2$thickness,
              mus = df2$lmus,
              thicknesses = c(65, 200, 500, 2000))

  # Run the LMP fit
  cat(n, "\n")
  fit = stan("/home/bbales2/CreepInfo/lumped.stan", data = data, warmup = 1000, iter = 1100, cores = 4)
  
  fits[[n]] = extract(fit, c("a", "b"))
}

list(n = map(fits, ~ .$a) %>% unlist,
     p = map(fits, ~ .$b) %>% unlist) %>%
  as.tibble %>%
  gather(name, value) %>%
  group_by(name) %>%
  summarize(m = mean(value),
            sd = sd(value))

  ggplot(aes(value)) +
  geom_histogram() +
  facet_grid(. ~ name)
