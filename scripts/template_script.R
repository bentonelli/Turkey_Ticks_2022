### Script to get started using mechanistic models ####


### Install packages ####
#You'll need to make sure that you install these packages before running the code
#Some need to be installed directly from the source

library(devtools)
install_github("eco4cast/EFIstandards")

library(lubridate)
library(ggplot2)
library(dplyr)
library(MMWRweek)
library(neon4cast)
library(xml2)
library(EFIstandards)
library(readr)


### Read in data,process ####

#Get data on sites
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ticks/master/Ticks_NEON_Field_Site_Metadata_20210928.csv")

#Reduce to metrics of interest
site_data_filtered <- select(site_data,c("field_site_id","field_site_county",
                                         "field_latitude","field_mean_annual_precipitation_mm",
                                         "field_mean_annual_temperature_C","field_avg_grean_increase_doy"))

#Get historical tick data
ticks_target <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)

### Analysis, visualization ####

#Look at tick data - some examples of how to look at things!
summary(ticks_target)
hist(year(ticks_target$time))
hist(log(ticks_target$amblyomma_americanum))
hist(ticks_target$mmwrWeek)
table(ticks_target$siteID)

#Create plot shoowing tick abundance at each site
cols_for_plot <- rainbow(length(unique(ticks_target$siteID)))
plot(ticks_target$mmwrWeek,ticks_target$amblyomma_americanum,cex=.8,pch=19,
     col=alpha(cols_for_plot[as.numeric(as.factor(ticks_target$siteID))],.8),xlab="MMWR Week",ylab="Tick Abundance")

#Gets the average across all sites, plot
tick_by_week <- ticks_target %>%
  group_by(mmwrWeek) %>%
  summarise(mean_tick_abund = mean(amblyomma_americanum))
points(tick_by_week,pch=15,col="black",type="l",lwd=4)

### Model creation ####

### Save predictions ####
### Submit forecast ####


