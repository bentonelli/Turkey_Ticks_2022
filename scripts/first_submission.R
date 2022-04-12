### Script to get data, make a basic plot ####
library(lubridate)
library(ggplot2)
library(dplyr)
library(MMWRweek)
library(neon4cast)
library(xml2)
library(devtools)
library(EFIstandards)
library(readr)
install_github("eco4cast/EFIstandards")
#Get data on sites
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ticks/master/Ticks_NEON_Field_Site_Metadata_20210928.csv")

#Reduce to metrics of interest
site_data_filtered <- select(site_data,c("field_site_id","field_site_county",
                                         "field_latitude","field_mean_annual_precipitation_mm",
                                         "field_mean_annual_temperature_C","field_avg_grean_increase_doy"))

#Get historical tick data
ticks_target <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)

#Look at tick data - basic 
summary(ticks_target)
hist(year(ticks_target$time))
hist(log(ticks_target$amblyomma_americanum))
hist(ticks_target$mmwrWeek)
table(ticks_target$siteID)
cols_for_plot <- rainbow(length(unique(ticks_target$siteID)))
plot(ticks_target$mmwrWeek,ticks_target$amblyomma_americanum,cex=.8,pch=19,
     col=alpha(cols_for_plot[as.numeric(as.factor(ticks_target$siteID))],.8),xlab="MMWR Week",ylab="Tick Abundance")

tick_by_week <- ticks_target %>%
  group_by(mmwrWeek) %>%
  summarise(mean_tick_abund = mean(amblyomma_americanum))

points(tick_by_week,pch=15,col="black",type="l",lwd=4)

#For each site, in each year, if data exists, fit loess function
all_site_pred <- list()
count <- 0
for (each_site in unique(ticks_target$siteID)){
  count <- count + 1
  site_pred <- c()
  for (each_year in 2014:2020){
    tick_yr_site <- filter(ticks_target,siteID == each_site & year(time)==each_year)
    if(nrow(tick_yr_site)>5){
      ls_pred <- loess(tick_yr_site$amblyomma_americanum~tick_yr_site$mmwrWeek)
      pred_ticks <- predict(ls_pred,newdata = 1:52) 
      site_pred <- rbind(site_pred,pred_ticks)
      plot(tick_yr_site$mmwrWeek,tick_yr_site$amblyomma_americanum)
      points(1:52,pred_ticks,type="l")
    }
  }
  all_site_pred[[count]] <- site_pred
}

#for each site, get mean, sd of each week
site_means <- c()
site_sds <- c()
predictions <- c()
for (i in 1:9){
  site_name <- unique(ticks_target$siteID)[i]
  site_in <- all_site_pred[[i]] 
  week_mean <- colMeans(site_in,na.rm = TRUE)
  #Change week means of < 0 to 1
  week_mean[is.na(week_mean) | week_mean <0] <- 1
  site_means <- rbind(week_mean,site_means)
  
  week_sd <- apply(site_in,2,sd,na.rm = TRUE)
  
  #Change site sds of NA to min sd for that row
  week_sd[is.na(week_sd)] <- quantile(week_sd,.25,na.rm = TRUE)
  
  week_97_5_upper <- week_mean + week_sd*1.96
  week_2_5_lower <- week_mean - week_sd*1.96
  
  week_2_5_lower[week_2_5_lower<0] <- 0
  
  site_sds <- rbind(site_sds,week_sd,week_97_5_upper,week_2_5_lower)
  
  week_add <- rep(1:52,4)
  add_site_pred <- as.data.frame(week_add)
  add_site_pred$siteID <- site_name
  
  add_site_pred$statistic <- NA
  add_site_pred$statistic[1:52] <- "mean"
  add_site_pred$statistic[53:104] <- "sd"
  add_site_pred$statistic[105:156] <- "Pred_interval_97.5"
  add_site_pred$statistic[157:208] <- "Pred_interval_02.5"
  
  add_site_pred$pred <- NA
  add_site_pred$pred[1:52] <- week_mean
  add_site_pred$pred[53:104] <- week_sd
  add_site_pred$pred[105:156] <- week_97_5_upper
  add_site_pred$pred[157:208] <- week_2_5_lower
  predictions <- rbind(predictions,add_site_pred)
}

weeks_needed <- c(10,11,12,13)
predictions_trim <- filter(predictions,week_add %in% weeks_needed)

predictions_trim$time <- MMWRweek::MMWRweek2Date(rep(2021,nrow(predictions_trim)),predictions_trim$week_add)
predictions_trim$forecast <- 1
predictions_trim$data_assimilation <- 0
predictions_trim$amblyomma_americanum <- predictions_trim$pred

predictions_trim <- select(predictions_trim,c("time","siteID","statistic",
                                              "forecast","data_assimilation",
                                              "amblyomma_americanum"))

#Make plots showing our predictions
for (each_site in unique(predictions_trim$siteID)){
  data_for_site <- filter(predictions_trim,siteID==each_site & statistic == "mean")
  sd_for_site <- filter(predictions_trim,siteID==each_site & statistic == "sd") 
  upper_for_site <- filter(predictions_trim,siteID==each_site & statistic == "Pred_interval_97.5") 
  lower_for_site <- filter(predictions_trim,siteID==each_site & statistic == "Pred_interval_02.5") 
  
  plot(week(data_for_site$time),data_for_site$amblyomma_americanum,
       main=each_site,type="l",lwd=5,col="dodgerblue4",ylim=c(0,max(upper)))
  
  polygon(c(min(week(data_for_site$time)):max(week(data_for_site$time)), 
            rev(c(min(week(data_for_site$time)):max(week(data_for_site$time))))),                # X-Coordinates of polygon
          c(upper_for_site$amblyomma_americanum, lower_for_site$amblyomma_americanum),  
          border = NA,# Y-Coordinates of polygon
          col = alpha("dodgerblue",.1)) 
  # Then plot the historical data
  ticks_target_site <- filter(ticks_target,siteID==each_site & mmwrWeek %in% unique(week(predictions_trim$time)))
  points(ticks_target_site$mmwrWeek,ticks_target_site$amblyomma_americanum,pch=19)
}

#Write to csv

write_csv(predictions_trim,"ticks-2021-03-07-UCLA_2022.csv")

#Create metadata file

forecast_file <- "ticks-2021-03-07-UCLA_2022.csv"

model_metadata = list(
  forecast = list(
    model_description = list(
      forecast_model_id =  "ticks", 
      name = "Combined loess model", 
      type = "empirical",  
      repository = "https://github.com/bentonelli/Turkey_Ticks_2022/scripts/template_script.R" 
    ),
    initial_conditions = list(
      status = "absent"
    ),
    drivers = list(
      status = "absent"
    ),
    parameters = list(
      status = "absent"
    ),
    random_effects = list(
      status = "absent"
    ),
    process_error = list(
      status = "data_driven",
      complexity = 2
    ),
    obs_error = list(
      status = "absent"
    )
  )
)
team_list <- list(list(individualName = list(givenName = "Ben", 
                                             surName = "Tonelli"),
                       organizationName = "UCLA",
                       electronicMailAddress = "btonelli@ucla.edu"),
                  list(individualName = list(givenName = "Graham", 
                                             surName = "Montgomery"),
                       organizationName = "UCLA"))

neon4cast::generate_metadata(forecast_file, team_list, model_metadata)

neon4cast::submit(forecast_file = "ticks-2021-03-07-UCLA_2022.csv",
                  metadata="ticks-2021-03-07-UCLA_2022.xml")

neon4cast::check_submission("ticks-2021-03-07-UCLA_2022.csv")
