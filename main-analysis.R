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

m1 <- brm(microcystin | cens(censored) ~ lake + t2(DOY, cYear, by = lake),
          data = mc, family = Gamma(link = "log"),
          warmup = 1000, iter = 3000, chains = 4, cores = CORES,
          control = list(adapt_delta = 0.95))
ms1 <- marginal_smooths(m1)

m2 <- brm(microcystin | cens(censored) ~ lake + t2(DOY, cYear, by = lake, k = c(10,10)),
          data = mc, family = Gamma(link = "log"),
          warmup = 1000, iter = 3000, chains = 4, cores = CORES,
          control = list(adapt_delta = 0.99, max_treedepth = 50))
ms2 <- marginal_smooths(m2)

m3 <- brm(microcystin | cens(censored) ~ lake + t2(DOY, cYear, k = c(10,10)),
          data = mc, family = Gamma(link = "log"),
          warmup = 1000, iter = 3000, chains = 4, cores = CORES,
          control = list(adapt_delta = 0.99, max_treedepth = 50))
ms3 <- marginal_smooths(m3)

m4 <- brm(microcystin | cens(censored) ~ lake + t2(DOY, cYear) +
              s(DOY, by = lake, m = 1) + s(cYear, by = lake, m = 1),
          data = mc, family = Gamma(link = "log"),
          warmup = 1000, iter = 3000, chains = 4, cores = CORES,
          control = list(adapt_delta = 0.99, max_treedepth = 50))
ms4 <- marginal_smooths(m4)

m5 <- brm(microcystin | cens(censored) ~ olake + t2(DOY, cYear) +
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

fitResponse <- with(mc, microcystin)
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
    geom_point(subset(mc, microcystin >= 0.16),
               mapping = aes(x = DOY, y = microcystin, colour = year, group = year),
               inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~ lake, scales = "free_y") +
    theme(legend.position = "top", legend.key.width = unit(2, "cm")) +
    labs(y = expression(Microcystin ~ mu*g~L^{-1})) +
    scale_colour_viridis(option = 'plasma', begin = 0, end = 0.9)
p2

## posterior predictive checks
m4.cms <- data.frame(mc = c(vapply(as.data.frame(p4), mean, numeric(1)),
                            vapply(as.data.frame(p4), median, numeric(1)),
                            dfd[['microcystin']]),
                     dataset = rep(c('est (mean)','Posterior Median','Observed'),
                                   each = nrow(dfd)))
brks <- c(0, 0.16, seq(0.26, ceiling(max(m4.cms[['mc']])), by = 0.2))

## summarising by the median seems to give better posterior predictive checks
m4.cms <- subset(m4.cms, subset = dataset != 'est (mean)')

## want to do the table for every row in p4 so need a function
brks <- c(0, 0.16, seq(0.36, ceiling(max(p4)), by = 0.2))
bin_fun <- function(x, brks) {
    table(cut(x, brks))
}
result <- apply(p4, 1L, bin_fun, brks = brks)
out <- apply(result, 1L, quantile, probs = c(0.05, 0.5, 0.95))
## brks_obs <- c(0, 0.16, seq(0.36, ceiling(max(m4.cms[['mc']])), by = 0.2))
tab_df <- data.frame(obs = unname(as.vector(table(cut(dfd[['microcystin']], brks)))[1:191]),
                     upper  = unname(out[3, 1:191]),
                     median = unname(out[2, 1:191]),
                     lower  = unname(out[1, 1:191]),
                     brks   = brks[1:191])

pp_dist_plt3 <- ggplot(tab_df[1:25, ]) +
    geom_point(aes(x = brks, y = obs, colour = 'Observed'), alpha = 1, size = 2) +
    geom_linerange(aes(x = brks, ymin = lower, ymax = upper), size = 1, colour = '#33bbee', alpha = 0.5) +
    geom_point(aes(x = brks, y = median, colour = 'Posterior'), size = 2, alpha = 0.5) +
    theme(legend.position = 'top') +
    scale_colour_manual(name = '', values = c('Observed' = '#cc3311', 'Posterior' = '#33bbee')) +
    labs(x = expression(Microcystin ~ mu*g~L^{-1}),
         y = "Frequency")
pp_dist_plt3

## Histogram posterior predictive check -----------------------------------------
pp_dist_plt <- ggplot(m4.cms, aes(x = mc, colour = dataset)) +
    geom_histogram(data = subset(m4.cms, subset = dataset == 'Posterior Median'),
                   mapping = aes(x = mc), breaks = brks, closed = 'left',
                   colour = '#ffffff', fill = '#33bbee') +
    geom_line(data = subset(m4.cms, subset = dataset == 'Observed'),
                  stat = 'bin', closed = 'right', breaks = brks,## size = 1,
                  show.legend = FALSE, colour = '#cc3311') +
    geom_point(data = subset(m4.cms, subset = dataset == 'Observed'),
               stat = 'bin', closed = 'right', breaks = brks, size = 0.75,
               show.legend = FALSE, colour = '#cc3311') +
    coord_cartesian(xlim = c(0,5)) +
    theme(legend.position = 'top') +
    labs(x = expression(Microcystin ~ mu*g~L^{-1}),
         y = "Frequency")
pp_dist_plt
## end --------------------------------------------------------------------------

## Frequency polygon of posterior predictive checks -----------------------------
pp_dist_plt2 <- ggplot(m4.cms, aes(x = mc, colour = dataset)) +
    geom_line(data = (m4.cms),
                  stat = 'bin', closed = 'right', breaks = brks) +
    geom_point(data = subset(m4.cms, subset = dataset == 'Posterior Median'),
               stat = 'bin', closed = 'right', breaks = brks, size = 0.75,
               show.legend = FALSE, colour = '#33bbee') +
    scale_colour_manual(name = "", values = c('#cc3311', '#33bbee')) +
    coord_cartesian(xlim = c(0,5)) +
    theme(legend.position = 'top') +
    labs(x = expression(Microcystin ~ mu*g~L^{-1}),
         y = "Frequency")
pp_dist_plt2
## end --------------------------------------------------------------------------

ggsave('posterior-predictive-check-freqpoly.pdf')

## posterior predictive checks on the number of censored values and number ------
## of values predicted values that exceed the limits.
## Compare against the observed metrics to see how wee we recreate the data
##  using the entire posterior
##
## Number of censored values first
cens_pp <- data.frame(prop_cens = c(rowMeans(p4 < 0.16)))
cens_obs <- mean(dfd[['microcystin']] < 0.16)

pp_cens_plt <- ggplot(cens_pp, aes(x = prop_cens)) +
    geom_histogram(bins = 30, colour = '#33bbee', fill = '#33bbee') +
    geom_vline(xintercept = cens_obs, colour = '#cc3311') +
    labs(x = "Proportion", y = "Frequency")

## next number of exceedances
levs <- c('> 0.3', '> 1.0', '> 1.6', '> 10', '> 20')
health_pp <- data.frame(prop = c(## rowMeans(p4 <  0.16),
                                 rowMeans(p4 >  0.3),
                                 rowMeans(p4 >  1.0),
                                 rowMeans(p4 >  1.6),
                                 rowMeans(p4 > 10.0),
                                 rowMeans(p4 > 20.0)),
                        threshold = ordered(rep(levs, each = nrow(p4)),
                                            levels = levs))

health_obs <- with(dfd, data.frame(prop = c(## mean(microcystin <   0.16),
                                            mean(microcystin >  0.3),
                                            mean(microcystin >  1.0),
                                            mean(microcystin >  1.6),
                                            mean(microcystin > 10.0),
                                            mean(microcystin > 20.0)),
                                   threshold = ordered(levs, levels = levs)))

pp_health_plt <- ggplot(health_pp, aes(x = prop)) +
    geom_histogram(bins = 50, colour = '#33bbee', fill = '#33bbee') +
    geom_vline(aes(xintercept = prop), data = health_obs, colour = '#cc3311') +
    facet_wrap(~ threshold, ncol = 3) +
    labs(x = "Proportion", y = "Frequency")

## put both of these plots together
plots <- align_plots(pp_dist_plt3, pp_health_plt, align = 'v', axis = 'l')
row1 <- plot_grid(plots[[1L]], pp_cens_plt, labels = 'auto', align = 'h', axis = 'ltb')
pp_plt <- plot_grid(row1, plots[[2L]], labels = c('','c'), rel_heights = c(1,1),
                    ncol = 1)
pp_plt
## end --------------------------------------------------------------------------

ggsave('posterior-predictive-checks.pdf', pp_plt, width = 7, height = 6.5)

## observed vs fitted values per lake -------------------------------------------
fitVals <- predict(m4, summary = TRUE, robust = TRUE)
fitVals <- transform(dfd, fitted = fitVals[,'Estimate'],
                     lake = factor(lake, levels = c('B','L','WW','P','K','C','D')))
fitVals <- droplevels(subset(fitVals, subset = lake != "D"))
siteVec <- c('Buffalo Pound','Last Mountain','Wascana','Pasqua','Katepwa','Crooked')
names(siteVec) <- levels(fitVals[['lake']])
quapelle_labeller <- as_labeller(siteVec)

fit_vs_obs_plt <- ggplot(fitVals, aes(x = fitted, y = microcystin,
                                      colour = as.logical(abs(censored)))) +
    geom_point(alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1, col = 'black', alpha = 0.3) +
    facet_wrap( ~ lake, labeller = quapelle_labeller, scales = 'free', ncol = 3) +
    theme(legend.position = 'top', legend.key.width = unit(1.5, 'cm')) +
    ## scale_colour_viridis_c(name = "Year", option = 'C') +
    labs(x = expression(Fitted   ~ Microcystin ~ (mu*g~L^{-1})),
         y = expression(Observed ~ Microcystin ~ (mu*g~L^{-1})),
         colour = NULL) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c("black", "red"), guid = "none") +
    theme(strip.background = element_rect(fill = "white")) ##  +
    ### coord_cartesian(xlim = c(0, 5), ylim = c(0,5))
fit_vs_obs_plt
## end --------------------------------------------------------------------------

## Plot of lake specific effects using variance components ----------------------
## lake specific effects
m4.summ <- summary(m4, waic = FALSE)$splines
lake_effs <- m4.summ[, c(1,3,4)]
lake_effs <- lake_effs[-(1:3), ]
lake_effs <- setNames(as.data.frame(lake_effs),
                      c('Estimate','Lower','Upper'))
lake_effs <- transform(lake_effs,
                       Lake = factor(rep(c('B','C','K','L','P','WW'), times = 2),
                                     labels = c('Buffalo Pound','Crooked','Katepwa',
                                                'Last Mountain','Pasqua','Wascana')),
                       Term = rep(c('Day of Year', 'Year'), each = 6))
lake_effs <- transform(lake_effs,
                       Lake = factor(Lake, levels = rev(c('Buffalo Pound','Last Mountain',
                                                          'Wascana','Pasqua','Katepwa','Crooked'))))

lake_specific_effs_plt <- ggplot(lake_effs, aes(x = Estimate, y = Lake)) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0) +
    geom_point() +
    facet_wrap(~ Term) +
    labs(y = NULL, x = 'Estimated standard deviation')
lake_specific_effs_plt
## end --------------------------------------------------------------------------

newd <- expand.grid(DOY = 120:243,
                    year = 2006:2016,
                    lake = c("B", "L", "WW", "P", "K", "C"))
newd <- transform(newd, cYear = year - refYear)

set.seed(1231323423)
pred4 <- posterior_predict(m4, newdata = newd)
saveRDS(pred4, "m4-posterior-predictions.rds")

set.seed(4587645786)
newd4 <- transform(newd,
                   mean = colMeans(pred4),
                   median = apply(pred4, 2L, quantile, probs = 0.5),
                   lower    = apply(pred4, 2L, quantile, probs = 0.025),
                   upper    = apply(pred4, 2L, quantile, probs = 0.975))
saveRDS(newd4, "m4-posterior-summaries.rds")


## Plot posterior predictions ---------------------------------------------------
post_pred_plt <- ggplot(newd4, aes(x = DOY, y = mean, colour = year, group = year)) +
    geom_line() +
    facet_wrap(~ lake, scales = 'free_y') + theme(legend.position = "top") +
    labs(y = expression(Microcystin ~ mu*g~L^{-1})) +
    scale_colour_viridis(option = 'plasma', begin = 0, end = 0.9)
post_pred_plt
## end --------------------------------------------------------------------------

## Plot posterior predictions and data ------------------------------------------
post_pred_data_plt <-
    ggplot(newd4, aes(x = DOY, y = mean, colour = year, group = year)) +
    geom_point(subset(dfd, microcystin >= 0.16),
               mapping = aes(x = DOY, y = microcystin, colour = year, group = year),
               inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~ lake, scales = "free_y") +
    theme(legend.position = "top", legend.key.width = unit(2, "cm")) +
    labs(y = expression(Microcystin ~ mu*g~L^{-1})) +
    scale_colour_viridis(option = 'plasma', begin = 0, end = 0.9)
post_pred_data_plt
## end --------------------------------------------------------------------------

## Plot of posterior probability a lake exceeds limits --------------------------
## How many days above x limit?
tmp4 <- transform(newd4,
                   `0.16` = colMeans(pred4 >= 0.16),
                   `1`    = colMeans(pred4 >= 1),
                   `5`    = colMeans(pred4 >= 5),
                  `10`    = colMeans(pred4 >= 10))
names(tmp4)[9:12] <- c(0.16, 1, 5, 10)
newd4t <- gather(tmp4, key = "limit", value = "prob", `0.16`, `1`, `5`, `10`)
newd4t <- dplyr::mutate(newd4t, limit = ordered(limit, levels = c(0.16, 1, 5, 10)))

prob_exceed_plt <- ggplot(newd4t, aes(x = DOY, y = prob, colour = year, group = year)) +
    geom_line() +
    facet_grid(limit ~ lake) +
    theme(legend.position = "top", legend.key.width = unit(2, "cm")) +
    scale_colour_viridis(option = 'plasma', begin = 0, end = 0.9)
prob_exceed_plt
## end --------------------------------------------------------------------------

## plot posterior mean fitted values vs observed MC -----------------------------
fitnew <- fitted(m4, newdata = newd, probs = c(0.025, 0.5, 0.975))
fitnew <- cbind(newd, as.data.frame(fitnew))

post_mean_fitted_plt <-
    ggplot(fitnew, aes(x = DOY, y = Estimate, colour = year, group = year)) +
    geom_point(subset(dfd, microcystin >= 0.16),
               mapping = aes(x = DOY, y = microcystin, colour = year, group = year),
               inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~ lake, scales = "free_y") +
    theme(legend.position = "top", legend.key.width = unit(2, "cm")) +
    scale_colour_viridis(option = 'plasma', begin = 0, end = 0.9) +
    labs(y = expression(Posterior ~ mean ~ microcystin ~ concentration ~ bgroup("[", mu*g~L^{-1}, "]")),
         x = "Day of year", colour = "")
post_mean_fitted_plt
## --- end ----------------------------------------------------------------------
