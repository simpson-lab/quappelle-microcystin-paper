## Script to analyse microcystin time series in Qu'Appelle lakes

## Load packages####
library("brms")
library("ggplot2")
library("dplyr")
library("bayesplot")
library("tidyr")
library('viridis')

## set the theme for plots
theme_set(theme_bw())

## load data
mc <- read.csv("./data/quappelle-lake-microcycstin-time-series.csv")
mc <- transform(mc,
                date = as.Date(date, "%Y-%m-%d"),
                censored = ifelse(microcystin < 0.16, -1, 0))
mc <- transform(mc,
                year = as.numeric(format(date,'%Y')),
                DOY = as.numeric(format(date,'%j')))

mc <- subset(mc, select = c(lake, year, DOY, date, microcystin, censored))
refYear <- with(mc, median(year))
mc <- transform(mc, cYear = year - refYear)

## order the data
mc <- mc[with(mc, order(date)),]
mc <- subset(mc, lake != "D")
mc <- transform(mc, olake = ordered(lake))

## Bayesian censored lognormal model fitted with Stan
## Refit?
CORES <- 4

m1 <- brm(micro_censored | cens(censored) ~ lake + t2(DOY, cYear, by = lake),
          data = mc, family = Gamma(link = "log"),
          warmup = 1000, iter = 3000, chains = 4, cores = CORES,
          control = list(adapt_delta = 0.95))
ms1 <- marginal_smooths(m1)

m2 <- brm(micro_censored | cens(censored) ~ lake + t2(DOY, cYear, by = lake, k = c(10,10)),
          data = mc, family = Gamma(link = "log"),
          warmup = 1000, iter = 3000, chains = 4, cores = CORES,
          control = list(adapt_delta = 0.99, max_treedepth = 50))
ms2 <- marginal_smooths(m2)

m3 <- brm(micro_censored | cens(censored) ~ lake + t2(DOY, cYear, k = c(10,10)),
          data = mc, family = Gamma(link = "log"),
          warmup = 1000, iter = 3000, chains = 4, cores = CORES,
          control = list(adapt_delta = 0.99, max_treedepth = 50))
ms3 <- marginal_smooths(m3)

m4 <- brm(micro_censored | cens(censored) ~ lake + t2(DOY, cYear) +
              s(DOY, by = lake, m = 1) + s(cYear, by = lake, m = 1),
          data = mc, family = Gamma(link = "log"),
          warmup = 1000, iter = 3000, chains = 4, cores = CORES,
          control = list(adapt_delta = 0.99, max_treedepth = 50))
ms4 <- marginal_smooths(m4)

m5 <- brm(micro_censored | cens(censored) ~ olake + t2(DOY, cYear) +
              s(DOY, by = olake) + s(cYear, by = olake),
          data = mc, family = Gamma(link = "log"),
          warmup = 1000, iter = 3000, chains = 4, cores = CORES,
          control = list(adapt_delta = 0.99, max_treedepth = 50))
ms5 <- marginal_smooths(m5)

summary(m1, waic = FALSE)

summary(m2, waic = FALSE)

summary(m3, waic = FALSE)

summary(m4, waic = FALSE)

summary(m5, waic = FALSE)

plot(ms1, stype = "raster", theme = theme_bw() + theme(legend.position = "top"))

plot(ms2, stype = "raster", theme = theme_bw() + theme(legend.position = "top"))

plot(ms3, stype = "raster", theme = theme_bw() + theme(legend.position = "top"))

plot(ms4, stype = "raster", theme = theme_bw() + theme(legend.position = "top"), ask = FALSE)

plot(ms5, stype = "raster", theme = theme_bw() + theme(legend.position = "top"), ask = FALSE)

WAIC(m1, m2, m3, m4, m5)

LOO(m1, m2, m3, m4, m5)

p1 <- posterior_predict(m1)
p2 <- posterior_predict(m2)
p3 <- posterior_predict(m3)
p4 <- posterior_predict(m4)
p5 <- posterior_predict(m5)

fitResponse <- with(mc, micro_censored)
fitResponse <- fitResponse[!is.na(fitResponse)]
take <- fitResponse >= 0.16
ppc_dens_overlay(fitResponse[take], p4[sample(nrow(p4), 50), take]) + theme_bw() + scale_x_log10()

newd <- expand.grid(DOY = 120:243,
                    year = 2006:2016,
                    lake = c("B", "L", "WW", "P", "K", "C"))
newd <- transform(newd, cYear = year - refYear)

pred4 <- posterior_predict(m4, newdata = newd)

newd4 <- transform(newd,
                   mean = colMeans(pred4),
                   median = apply(pred4, 2L, quantile, probs = 0.5),
                   lower    = apply(pred4, 2L, quantile, probs = 0.025),
                   upper    = apply(pred4, 2L, quantile, probs = 0.975))


## These aren't the *final* figures in the paper; quick versions here to
##  just visualise the trends
p1 <- ggplot(newd4, aes(x = DOY, y = mean, colour = year, group = year)) +
    geom_line() +
    facet_wrap(~ lake, scales = 'free_y') + theme(legend.position = "top") +
    labs(y = expression(Microcystin ~ mu*g~L^{-1})) +
    scale_colour_viridis(option = 'plasma', begin = 0, end = 0.9)
p1

p2 <- ggplot(newd4, aes(x = DOY, y = mean, colour = year, group = year)) +
    geom_point(subset(mc, micro_censored >= 0.16),
               mapping = aes(x = DOY, y = micro_censored, colour = year, group = year),
               inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~ lake, scales = "free_y") +
    theme(legend.position = "top", legend.key.width = unit(2, "cm")) +
    labs(y = expression(Microcystin ~ mu*g~L^{-1})) +
    scale_colour_viridis(option = 'plasma', begin = 0, end = 0.9)
p2

