## Load packages####
library("brms")
library("ggplot2")
library("bayesplot")
library("tidyr")
library('viridis')
library("readr")
library("dplyr")
library('cowplot')

## set the theme for plots
theme_set(theme_bw())

## load data
mc <- read_csv("./data/microcystin-data.csv")
bg <- read_csv("./data/other-variables.csv")

lake_levels <- c("B", "L", "WW", "P", "K", "C")

names(mc) <- c("lake", "date", "year", "DOY", "cYear", "microcystin")
mc <- mutate(mc,
             lake = factor(lake, levels = lake_levels),
             date = as.Date(date, "%Y-%m-%d"),
             censored = ifelse(microcystin < 0.16, -1, 0))
mc <- mutate(mc,
             year = as.numeric(format(date,'%Y')),
             DOY = as.numeric(format(date,'%j')))

mc <- select(mc, lake, year, DOY, date, microcystin, censored)
refYear <- with(mc, median(year))
mc <- mutate(mc, cYear = year - refYear)

bg <- mutate(bg, lake = factor(lake, level = lake_levels))

## order the data
mc <- mc[with(mc, order(date)),]
mc <- filter(mc, lake != "D")
mc <- mutate(mc, olake = ordered(lake))
## add a year factor
mc <- mutate(mc, fYear = factor(year))

##
df <- left_join(mc, select(bg, lake, year, DOY, temp_surface, TDN_ug_L, TDP_ug_L),
                by = c("lake", "year", "DOY"))
df <- mutate(df, log_tdn = log(TDN_ug_L), log_tdp = log(TDP_ug_L))

## Bayesian censored lognormal model fitted with Stan
#do a lake specific model to look model MC time series
CORES <- 4

REFIT <- FALSE
if (REFIT) {
    ## model with DOY, year, temp, and TDN, all lakes
    m1 <- brm(microcystin | cens(censored) ~ lake + t2(DOY, year, by = lake) +
                  s(temp_surface) + s(log_tdn),
              data = df, family = Gamma(link = "log"),
              cores = CORES, warmup = 1000, iter = 2000, chains = 4,
              control = list(adapt_delta = 0.999, max_treedepth = 25))
    
    ## model with DOY, year, temp, and TDN, all lakes
    m2 <- brm(microcystin | cens(censored) ~ lake +
                  t2(DOY, year, by = lake) +
                  s(temp_surface, by = lake) +
                  s(log_tdn, by = lake),
              data = df, family = Gamma(link = "log"),
              cores = CORES, warmup = 1000, iter = 2000, chains = 4,
              control = list(adapt_delta = 0.9999, max_treedepth = 25))
    
    ## Log TDP
    ## model with DOY, year, temp, and TDN, all lakes
    m3 <- brm(microcystin | cens(censored) ~ lake + t2(DOY, year, by = lake) +
                  s(temp_surface) + s(log_tdp),
              data = df, family = Gamma(link = "log"),
              cores = CORES, warmup = 1000, iter = 2000, chains = 4,
              control = list(adapt_delta = 0.999, max_treedepth = 25))
    
    ## model with DOY, year, temp, and TDN, all lakes
    m4 <- brm(microcystin | cens(censored) ~ lake +
                  t2(DOY, year, by = lake) +
                  s(temp_surface, by = lake) +
                  s(log_tdp, by = lake),
              data = df, family = Gamma(link = "log"),
              cores = CORES, warmup = 1000, iter = 2000, chains = 4,
              control = list(adapt_delta = 0.9999, max_treedepth = 25))
    
    saveRDS(m1, "m1-tdn.rds")
    saveRDS(m2, "m2-tdn.rds")
    saveRDS(m3, "m3-tdp.rds")
    saveRDS(m4, "m4-tdp.rds")
} else {
    m1 <- readRDS("m1-tdn.rds")
    m2 <- readRDS("m2-tdn.rds")
    m3 <- readRDS("m3-tdp.rds")
    m4 <- readRDS("m4-tdp.rds")
}


ms1 <- marginal_smooths(m1)

saveRDS(ms1, 'temp-tdn-m1-marginal-smooths.rds')

plot(ms1)

ms2 <- marginal_smooths(m2)

plot(ms2)

ms3 <- marginal_smooths(m3)

plot(ms3)

ms4 <- marginal_smooths(m4)

plot(ms4)

## m1 is with TDN Reviewer asked whether the seasonal smooth changed when
## we included the covariates in the model
ms1doy <- ms1[["mu: t2(DOY,year,by=lake)"]]

doy_plt_with_tdn <- ggplot(ms1doy, aes(x = DOY, y = estimate__, colour = year, group = year)) +
    ## geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.1) +
    geom_line() +
    facet_wrap(~ lake, scales = "free_y") +
    theme(legend.position = "top", legend.key.width = unit(2, "cm")) +
    scale_colour_viridis(option = 'plasma', begin = 0, end = 0.9) +
    labs(x = "Day of Year", y = "Centred effect", colour = "Year")
doy_plt_with_tdn

## load m4 and ms4 as those were the original models described in paper
m4 <- readRDS("../m4-te-with-doy-and-year-departure-smooths-default-k.rds")
ms4 <- readRDS("../ms4-te-with-doy-and-year-departure-smooths-default-k.rds")


## m1 is with TDN
ms1temp <- ms1[["mu: s(temp_surface)"]]
ms1tdn <- ms1[["mu: s(log_tdn)"]]

plt1temp <- ggplot(ms1temp, aes(x = temp_surface, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.1) +
    geom_line() +
    labs(x = expression(Surface ~ temperature ~ (degree*C)),
         y = "Centred effect")

plt1tdn <- ggplot(ms1tdn, aes(x = log_tdn, y = estimate__)) +
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.1) +
    geom_line() +
    labs(x = expression(log[e](TDN)),
         y = "Centred effect")

covar_temp_tdn_plts <- plot_grid(plt1temp, plt1tdn, ncol = 2,
                                 align = 'vh', axis = 'tb', labels = 'auto')
covar_temp_tdn_plts

ggsave('covariate-effects-temp-tdn.pdf', covar_temp_tdn_plts, width = 9, height = 4.5)

## Load m4 from the main paper:
main4 <- readRDS("../m4-te-with-doy-and-year-departure-smooths-default-k.rds")

## lake specific effects
main4.summ <- summary(main4, waic = FALSE)$splines
lake_effs <- main4.summ[, c(1,3,4)]
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
    labs(y = NULL, x = 'Estimated standard deviation') +
    theme(strip.background = element_rect(fill = "white"))
lake_specific_effs_plt

row1 <- lake_specific_effs_plt
# first align the top-row plot (row1) with the left-most plot of the
# bottom row (plt1temp)
plots <- align_plots(row1, plt1temp, align = 'v', axis = 'l')
row2 <- plot_grid(plots[[2]], plt1tdn, ncol = 2,
                  align = 'h', axis = 'tb', labels = c('B.','C.'))

figure_3 <- plot_grid(plots[[1]], row2, nrow = 2, labels = c('A.', '', ''),
                      rel_heights = c(1,1))
figure_3

ggsave('./figure_3.pdf', width = 8, height = 6)
