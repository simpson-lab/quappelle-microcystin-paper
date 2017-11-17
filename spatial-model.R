## Spatial smooth of Heather's data

## Load packages####
library("brms")
library("ggplot2")
library("bayesplot")
library("cowplot")
theme_set(theme_bw())
library("tidyr")
library("readr")
library('viridis')

## load data
survey <- read_csv("./data/sask-lakes-survey-microcystin-data.csv")
survey <- transform(survey,
                    censored = ifelse(MC < 0.16, -1, 0))

## Refit?
CORES <- 4

msos <- brm(bf(MC | cens(censored) ~ s(longitude, latitude, bs = "sos", k = 60)),
            data = survey,
            family = Gamma(link = "log"),
            warmup = 1500, iter = 5500, chains = 4, cores = CORES,
            control = list(adapt_delta = 0.999999, max_treedepth = 100),
            seed = 42)
msos.ms <- marginal_smooths(msos)

summary(msos, waic = FALSE)
plot(msos, stype = "raster")

mtp <- brm(bf(MC | cens(censored) ~ t2(longitude, latitude, k = 9)),
           data = survey,
           family = Gamma(link = "log"),
           warmup = 1500, iter = 5500, chains = 4, cores = CORES,
           control = list(adapt_delta = 0.99999, max_treedepth = 80),
           seed = 24)
mtp.ms <- marginal_smooths(mtp)

summary(mtp, waic = FALSE)
plot(mtp, stype = "raster")

WAIC(msos, mtp)

LOO(msos, mtp) ## reloo = TRUE ended up throwing an error - invalid 'times' argument
## results not too different to the WAIC, the SOS model is slightly better

pdata <- with(survey, expand.grid(latitude = seq(min(latitude), max(latitude), length = 100),
                                  longitude = seq(min(longitude), max(longitude), length = 100)))

fitvals  <- fitted(mtp,  newdata = pdata, scale = 'response')
predvals <- predict(mtp, newdata = pdata, scale = 'response')

fitdata <- cbind(pdata, fitvals)
preddata <- cbind(pdata, predvals)

basePlt <- ggplot(fitdata, aes(x = longitude, y = latitude)) + coord_quickmap() +
    scale_fill_viridis(option = 'plasma') +
    theme(legend.position = 'top', legend.key.width = unit(1, 'cm')) +
    labs(x = 'Longitude', y = 'Latitude')
basePredPlt <- ggplot(preddata, aes(x = longitude, y = latitude)) + coord_quickmap() +
    scale_fill_viridis(option = 'plasma') +
    theme(legend.position = 'top', legend.key.width = unit(1, 'cm')) +
    labs(x = 'Longitude', y = 'Latitude')
sampleLayer <- geom_point(data = survey, mapping = aes(x = longitude, y = latitude))

p1 <- basePlt + geom_raster(aes(fill = Estimate)) + geom_contour(aes(z = Estimate)) +
    sampleLayer
p1l <- basePlt + geom_raster(aes(fill = `2.5%ile`)) + sampleLayer
p1u <- basePlt + geom_raster(aes(fill = `97.5%ile`)) + sampleLayer
p1se <- basePlt + geom_raster(aes(fill = `Est.Error`)) + sampleLayer
p2 <- basePredPlt + geom_raster(aes(fill = Estimate)) + sampleLayer
p2l <- basePredPlt + geom_raster(aes(fill = `2.5%ile`)) + sampleLayer
p2u <- basePredPlt + geom_raster(aes(fill = `97.5%ile`)) + sampleLayer

pp <- plot_grid(p1, p2, p2l, p2u, ncol = 2, align = 'hv')
pp

vals <- c(0.01, 0.1, 0.5, 1, 5, 10, 20, 40, 50)
sp1 <- ggplot(survey, aes(x = longitude, y = latitude, colour = log10(MC))) +
    coord_quickmap() +
    scale_colour_viridis(name = expression(mu*g~L^{-1}),
                         option = 'plasma',
                         breaks = round(log10(vals), 2),
                         labels = vals,
                         end = 0.9, begin = 0.1) +
    theme(legend.position = 'top', legend.key.width = unit(1, 'cm')) +
    labs(x = 'Longitude', y = 'Latitude') +
    geom_point(size = 1.5)
sp1

sp2 <- ggplot(survey, aes(x = MC)) +
    geom_histogram(closed = 'left', breaks = c(0, 0.16, seq(1, 45, by = 1))) +
    labs(x = expression(Micocystin ~ (mu*g~L^{-1})), y = "Frequency")
sp2

plts <- plot_grid(sp1, sp2, ncol = 2, align = 'h', axis = 'b',
                  labels = 'auto')
plts

