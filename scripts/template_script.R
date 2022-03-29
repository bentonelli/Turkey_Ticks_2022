### Script to get data ####
library(lubridate)
library(ggplot2)
library(dplyr)

#Get data on sites
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ticks/master/Ticks_NEON_Field_Site_Metadata_20210928.csv")

#Get historical tick data
ticks_target <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)

#Look at tick data - basic 
summary(ticks_target)
hist(year(ticks_target$time))
hist(ticks_target$amblyomma_americanum)
hist(ticks_target$mmwrWeek)
table(ticks_target$siteID)
cols_for_plot <- rainbow(length(unique(ticks_target$siteID)))
plot(ticks_target$mmwrWeek,ticks_target$amblyomma_americanum,cex=.8,pch=19,
     col=alpha(cols_for_plot[as.numeric(as.factor(ticks_target$siteID))],.8),xlab="MMWR Week",ylab="Tick Abundance")

tick_by_week <- ticks_target %>%
  group_by(mmwrWeek) %>%
  summarise(mean_tick_abund = mean(amblyomma_americanum))

points(tick_by_week,pch=15,col="black",type="l",lwd=4)
