## GAM models for supplemental figures in Hayes et al L&O: Letters

## Nicole Hayes

## Package
library("mgcv")
library("lubridate")
library("rLakeAnalyzer")
library("plyr") # for revalue()

##Load profile data
profiles <- read.csv("./data/quappelle-lake-profile-data.csv")
profiles <- transform(profiles, Date= as.Date(as.character(Date)))
names(profiles)[names(profiles) == "Day_of_Year"] <- "DOY"

#Surface Temperature model; Fig. S3####
surf <- subset(profiles, Depth == "0")

mod <- gam(Temperature ~ Lake + te(Year, DOY, by = Lake), data = surf,
           method = "REML", select = TRUE, family = gaussian())

#Schmidt Stability; Fig. S4####
head(profiles)

###lake by lake....
#Katepwa####
Katepwa <- subset(profiles, Lake == "K")
Katepwa <- subset(Katepwa, Date != "2008-08-01") #No surface temperature
Katepwa <- subset(Katepwa, Depth != "0.5")

#Create formula for area curve for Katepwa, data from capacity curve pdf
elev_KK <- c(455.2,458,460,462,464,466,468,470,472,473,476,478,480,481.6)
area_KK <- c(0,496,704,816,930,1028,1120,1184,1280,1384,1480,1572,1770,1940)# area is in ha
area_KK <- area_KK*10000 #convert to m2

## make a single dataframe of all the data
lake <- rep(c("Katepwa"),
            times = c(14))
spline <- data.frame(lake, elev_KK, area_KK)
names(spline) <- c("lake", "elevation", "area")

##assign a surface area to each depth.
spl <- with(Katepwa, split(Katepwa, Date))

fn <- function(df) {
  #get the max depth measured for every sampling date
  count <- which.max(df[, "Depth"])
  depth <- with(df, Depth[count])
  #Now convert those depths to lake elevations
  base <- which.min(spline[, "elevation"])
  elev <- with(spline, elevation[base])
  elevation <- seq(from = (elev + 1),  to = (elev + 1) + depth, by = 1) ##assume the deepest sample we take is 1m from the bottom of the lake
  #create a spline describing the elevation to lake area realtionship
  spFun <- with(spline, splinefun(elevation, area, method = "natural"))
  ##now apply the spline to the elevations sampled for each date
  SA_sampleDepth <- spFun(elevation)
  with(df, data.frame(Lake = Lake, Date = Date, Year = Year, Depth = Depth,
                      Temperature = Temperature, SA = rev(SA_sampleDepth)))
}

out <- lapply(spl, fn)
Katepwa_SA <- do.call("rbind", out)

#Crooked####
Crooked <- subset(profiles, Lake == "C")
Crooked<- subset(Crooked, Date != "1995-05-15") # No profile data
Crooked<- subset(Crooked, Date != "2001-05-16") # No temperature data
Crooked<- subset(Crooked, Date != "2006-08-29") # Missing profile depths
Crooked<- subset(Crooked, Date != "2007-01-20") # Remove; winter sampling date
Crooked<- subset(Crooked, Date != "2007-02-20") # Remove; winter sampling date
Crooked<- subset(Crooked, Date != "2011-05-11") #Missing surface sample
Crooked<- subset(Crooked, Depth != "0.5") # Remove half meter sampling depths
Crooked<- subset(Crooked, Depth != "1.5")
Crooked<- subset(Crooked, Depth != "2.5")
Crooked<- subset(Crooked, Depth != "3.5")
Crooked<- subset(Crooked, Depth != "4.5")
Crooked<- subset(Crooked, Depth != "5.5")
Crooked<- subset(Crooked, Depth != "6.5")
Crooked<- subset(Crooked, Depth != "7.5")
Crooked<- subset(Crooked, Depth != "8.5")
Crooked<- subset(Crooked, Depth != "9.5")
Crooked<- subset(Crooked, Depth != "10.5")
Crooked <- Crooked[complete.cases(Crooked[, 3]),]

#Create formula for area curve for Katepwa, data from capacity curve pdf
elev_C <- c(435, 436, 438, 440, 442, 444, 446, 448, 450, 452, 454)
area_C <- c(0, 250, 355, 523, 648, 762, 870, 980, 1185, 1620, 2225)# area is in ha
area_C <- area_C*10000 #convert to m2

## make a single dataframe of all the data
lake <- rep(c("Crooked"),
            times = c(11))
spline <-data.frame(lake, elev_C, area_C)
names(spline) <- c("lake", "elevation", "area")

##assign a surface area to each depth....
spl <- with(Crooked, split(Crooked, Date))

fn <- function(df) {
  #get the max depth measured for every sampling date
  count <- which.max(df[, "Depth"])
  depth <- with(df, Depth[count])
  #Now convert those depths to lake elevations
  base <- which.min(spline[, "elevation"])
  elev <- with(spline, elevation[base])
  elevation <- seq(from = (elev + 1),  to = (elev + 1) + depth, by = 1) ##assume the deepest sample we take is 1m from the bottom of the lake
  #create a spline describing the elevation to lake area realtionship
  spFun <- with(spline, splinefun(elevation, area, method = "natural"))
  ##now apply the spline to the elevations sampled for each date
  SA_sampleDepth <- spFun(elevation)
  with(df, data.frame(Lake = Lake, Date = Date, Year = Year, Depth = Depth,
                      Temperature = Temperature, SA = rev(SA_sampleDepth)))
}

out <- lapply(spl, fn)
C_SA <- do.call("rbind", out)

#Pasqua####
Pasqua <- subset(profiles, Lake == "P")
Pasqua <- subset(Pasqua, Date != "2009-05-26") # No surface sample
Pasqua <- subset(Pasqua, Date != "2006-07-05") # Missing profile depths
Pasqua <- subset(Pasqua, Depth != "0.5")
Pasqua <- subset(Pasqua, Depth != "1.5")
Pasqua <- subset(Pasqua, Depth != "2.5")
Pasqua <- subset(Pasqua, Depth != "3.5")
Pasqua <- subset(Pasqua, Depth != "4.5")
Pasqua <- subset(Pasqua, Depth != "5.5")
Pasqua <- subset(Pasqua, Depth != "6.5")
Pasqua <- subset(Pasqua, Depth != "7.5")
Pasqua <- subset(Pasqua, Depth != "8.5")
Pasqua <- subset(Pasqua, Depth != "9.5")
Pasqua <- subset(Pasqua, Depth != "10.5")
Pasqua <- subset(Pasqua, Depth != "11.5")
Pasqua <- subset(Pasqua, Depth != "12.5")
Pasqua <- subset(Pasqua, Depth != "17.5")

#Create formula for area curve for Pasqua, data from capacity curve pdf
elev <- c(463,465,467,469,471,472,473,474,475,476,477,478,480,481.6)
area <- c(0,175,428,618,745,794,837,875,922,1032,1260,1585,2260,2570)# area is in ha
area <- area*10000 #convert to m2

## make a single dataframe of all the data
spline <-data.frame(elev, area)
names(spline) <- c("elevation", "area")

##assign a surface area to each depth.
spl <- with(Pasqua, split(Pasqua, Date))

fn <- function(df) {
  #get the max depth measured for every sampling date
  count <- which.max(df[, "Depth"])
  depth <- with(df, Depth[count])
  #Now convert those depths to lake elevations
  base <- which.min(spline[, "elevation"])
  elev <- with(spline, elevation[base])
  elevation <- seq(from = (elev + 1),  to = (elev + 1) + depth, by = 1) ##assume the deepest sample we take is 1m from the bottom of the lake
  #create a spline describing the elevation to lake area realtionship
  spFun <- with(spline, splinefun(elevation, area, method = "natural"))
  ##now apply the spline to the elevations sampled for each date
  SA_sampleDepth <- spFun(elevation)
  with(df, data.frame(Lake = Lake, Date = Date, Year = Year, Depth = Depth,
                      Temperature = Temperature, SA = rev(SA_sampleDepth)))
}

out <- lapply(spl, fn)
Pasqua_SA <- do.call("rbind", out)

#Now onto Buffalo Pound####
BP <- subset(profiles, Lake == "B")
BP <- subset(BP, Date != "2006-07-17") # No surface sample
BP <- subset(BP, Depth != "0.5")
BP <- subset(BP, Depth != "1.5")
BP <- subset(BP, Depth != "2.5")
BP <- subset(BP, Depth != "2.53")
BP <- subset(BP, Depth != "3.5")
BP <- subset(BP, Depth != "4.5")

#Create formula for area curve, data from capacity curve pdf
elev <-c(503.8,505,506,507,508,509,510,511,512,513,513.6)
area <-c(0,785,1512,1963,2250,2855,3030,3200,3370,3475,3680)
area <- area*10000 # convert to m2

## make a single dataframe of all the data
spline <-data.frame(elev, area)
names(spline) <- c("elevation", "area")

##assign a surface area to each depth....
spl <- with(BP, split(BP, Date))

fn <- function(df) {
  #get the max depth measured for every sampling date
  count <- which.max(df[, "Depth"])
  depth <- with(df, Depth[count])
  #Now convert those depths to lake elevations
  base <- which.min(spline[, "elevation"])
  elev <- with(spline, elevation[base])
  elevation <- seq(from = (elev + 1),  to = (elev + 1) + depth, by = 1) ##assume the deepest sample we take is 1m from the bottom of the lake
  #create a spline describing the elevation to lake area realtionship
  spFun <- with(spline, splinefun(elevation, area, method = "natural"))
  ##now apply the spline to the elevations sampled for each date
  SA_sampleDepth <- spFun(elevation)
  with(df, data.frame(Lake = Lake, Date = Date, Year = Year, Depth = Depth,
                      Temperature = Temperature, SA = rev(SA_sampleDepth)))
}

out <- lapply(spl, fn)
BP_SA <- do.call("rbind", out)

#Now onto Last Mountain####
LM <- subset(profiles, Lake == "L")
LM <- subset(LM, Date != "1995-05-12") # Only surface sample
LM <- subset(LM, Date != "1996-05-08") # Only samples to 2-m
LM <- subset(LM, Date != "1997-07-04") # Only surface sample
LM <- subset(LM, Date != "2015-05-29") # No profile data
LM <- LM[complete.cases(LM[, 6]),]

#Create formula for area curve, data from capacity curve pdf
elev_LM <- c(458.7,462,466,470,474,476,478,480,482,484,486,488,490,492,493.8)
area_LM <- c(0,40,370,1000,2150,3010,4020,5550,7750,10600,13550,16080,18100,23100,28530)
area_LM <- area_LM*10000

## make a single dataframe of all the data
spline <-data.frame(elev, area)
names(spline) <- c("elevation", "area")

##assign a surface area to each depth....
spl <- with(LM, split(LM, Date))

fn <- function(df) {
  #get the max depth measured for every sampling date
  count <- which.max(df[, "Depth"])
  depth <- with(df, Depth[count])
  #Now convert those depths to lake elevations
  base <- which.min(spline[, "elevation"])
  elev <- with(spline, elevation[base])
  elevation <- seq(from = (elev + 1),  to = (elev + 1) + depth, by = 1) ##assume the deepest sample we take is 1m from the bottom of the lake
  #create a spline describing the elevation to lake area realtionship
  spFun <- with(spline, splinefun(elevation, area, method = "natural"))
  ##now apply the spline to the elevations sampled for each date
  SA_sampleDepth <- spFun(elevation)
  with(df, data.frame(Lake = Lake, Date = Date, Year = Year, Depth = Depth,
                      Temperature = Temperature, SA = rev(SA_sampleDepth)))
}

out <- lapply(spl, fn)
LM_SA <- do.call("rbind", out)
LM_SA <- LM_SA[complete.cases(LM_SA[, 5]),]

#Now onto Wascana####
WW <- subset(profiles, Lake == "WW")
WW <- subset(WW, Depth != "0.5")
WW <- subset(WW, Depth != "1.5")
WW <- subset(WW, Depth != "2.5")
WW <- subset(WW, Depth != "3.5")
WW <- subset(WW, Depth != "4.5")

#Create formula for area curve, data from capacity curve pdf
elev <- c(564, 565,565.2,566, 567, 568, 569, 570, 571)
area <- c(0, 2.4697,25.0590, 30.7538, 34.9830, 38.5631, 46.8523, 116.8272, 206.4147)
area <- area *10000

## make a single dataframe of all the data
spline <-data.frame(elev, area)
names(spline) <- c("elevation", "area")

##assign a surface area to each depth....
spl <- with(WW, split(WW, Date))

fn <- function(df) {
  #get the max depth measured for every sampling date
  count <- which.max(df[, "Depth"])
  depth <- with(df, Depth[count])
  #Now convert those depths to lake elevations
  base <- which.min(spline[, "elevation"])
  elev <- with(spline, elevation[base])
  elevation <- seq(from = (elev + 1),  to = (elev + 1) + depth, by = 1) ##assume the deepest sample we take is 1m from the bottom of the lake
  #create a spline describing the elevation to lake area realtionship
  spFun <- with(spline, splinefun(elevation, area, method = "natural"))
  ##now apply the spline to the elevations sampled for each date
  SA_sampleDepth <- spFun(elevation)
  with(df, data.frame(Lake = Lake, Date = Date, Year = Year, Depth = Depth,
                      Temperature = Temperature, SA = rev(SA_sampleDepth)))
}

out <- lapply(spl, fn)
WW_SA <- do.call("rbind", out)

#Join up all the dataframes####
data <- rbind(Katepwa_SA, Pasqua_SA, BP_SA, WW_SA, LM_SA, C_SA)
data <- transform(data, ID = paste(Lake, Date, sep = "_"))

#feed it into LakeAnalyzer and calc Schmidt
#Split, Apply, Combine with a function to calculate buoyancy####
###
spl <- split(data, data$ID)

fn <- function(df) {
  out <- if (nrow(df) > 0) {
    wtr = with(df, Temperature)
    depths = with(df, Depth)
    bthA = with(df, SA)
    bthD = with(df, Depth)
    SS <- schmidt.stability(wtr, depths, bthA, bthD, sal = wtr * 0)
    with(df, data.frame(ID = ID[1], Lake = Lake[1], Date = Date[1],  SS = SS))
  } else {
    data.frame(ID = NA, Lake = NA, Date = NA,  depths = NA, wtr = NA, buoyancy = NA)
  }
  out
}

out <- lapply(spl, fn)
out <- do.call("rbind", out)
out <- transform(out, Year= as.numeric(format(Date, format = '%Y')))
out <- transform(out, DOY= yday(Date))
out <- subset(out, DOY >= 120 & DOY <= 243)
out <- transform(out, Name = revalue(Lake, c("B"="Buffalo Pound", "L" = "Last Mountain", "WW"= "Wascana",
                                             "P"= "Pasqua", "K"="Katepwa", "C" = "Crooked")))

#GAM for Schmidt Stability
modSS <- gam(SS ~ Lake + te(DOY, Year, by = Lake), data = out,
             method = "REML", select = TRUE, family = gaussian(), na.action = na.exclude)

summary(modSS)

#Canthaxanthin surface model; Fig. S5 ####
df <- read.csv("./data/quappelle-lake-time-series-data.csv")
df <- transform(df, Lake = factor(Lake, levels = c("B", "L", "WW", "P","K", "C")))

modCanth <- gam(Surface.canthaxanthin  ~ Lake + ti(Day.of.Year) + ti(Year) +
                    ti(Day.of.Year, Year, by = Lake) + ti(Day.of.Year, Year) +
                    ti(Day.of.Year, by = Lake) + ti(Year, by = Lake), data = df,
                method = "REML", select = TRUE, family = tw(), na.action = na.exclude)

summary(modCanth)

#Total dissolved nitrogen; Fig. S6####
modTDN <- gam(Total.Dissolved.Nitrogen ~ Lake + ti(Day.of.Year) + ti(Year) +
                  ti(Day.of.Year, Year, by = Lake) + ti(Day.of.Year, Year) +
                  ti(Day.of.Year, by = Lake) + ti(Year, by = Lake), data = df,
              method = "REML", select = TRUE, family = Gamma())

summary(modTDN)

#Total dissolved phosphorus; Fig. S7####
modTDP <- gam(Total.Dissolved.Phosphorus  ~ Lake + ti(Day.of.Year) + ti(Year) +
                  ti(Day.of.Year, Year, by = Lake) + ti(Day.of.Year, Year) +
                  ti(Day.of.Year, by = Lake) + ti(Year, by = Lake), data = df,
              method = "REML", select = TRUE, family = Gamma())

summary(modTDP)
