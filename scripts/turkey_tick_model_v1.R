# Function to run tick population model, V1.0 
library(dplyr)
library(lubridate)

#Read in our training data
ticks_target <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)
ticks_target <- filter(ticks_target,siteID == "SCBI")
#Read in our environmental data
temp_data <- read.csv("data/temp_data/Weekly Average Temperature/ SCBI .csv")
temp_data$year <- as.numeric(substr(temp_data$Week, 1, 4)) 
temp_data$week_num <- as.numeric(substr(temp_data$Week, 6, 7)) 

temp_data_avg <- temp_data %>% 
  group_by(week=temp_data$week_num) %>%
  summarise(avg_week_temp = mean(tempMean))
temp_data_avg$stand_temp <- (temp_data_avg$avg_week_temp - mean(temp_data_avg$avg_week_temp))

precip_data <- read.csv("data/Combined_weekly_precip_data.csv")
precip_data <- filter(precip_data,siteID == "SCBI")
precip_data_avg <- precip_data %>%
  group_by(week) %>%
  summarise(mean_week_precip = mean(sum_weekly_total_precip_bySite))
precip_data_avg$stand_prec <- precip_data_avg$mean_week_precip - mean(precip_data_avg$mean_week_precip)
plot(precip_data_avg$week,precip_data_avg$stand_prec)

turkey_tick_model_v1 <- function(temp_data = temp_data_avg,precip_data = precip_data_avg,
                                 br_larvae_molt=.001,br_nymph_molt=.001,br_adult_breed=.0185,
                                 max_larvae_molt=.01,max_nymph_molt=.01,
                                 exp_temp_param_larv_molt=0.0000,exp_temp_param_nymph_molt = 0.0001,
                                 exp_precip_param_larv_molt=0.00001,exp_precip_param_nymph_molt = 0.00001,
                                 exp_temp_param_adult_breed=0.0000,
                                days=365,starting_nymphs = 10,starting_adults=10){
  
  # Setting up the # days the simulation runs for
  days <- 1:365
  #Starting population sizes
  eggs <- 0
  larvae <- 0 #by life stage
  larvae_molt <- 0
  nymphs <- starting_nymphs
  nymph_molt <- 0
  adults <- starting_adults
  
  #Setting up record keeping
  tick_record <- matrix(NA, nrow=6,ncol=length(days)+1)
  tick_record[1:6,1] <- c(eggs,larvae,larvae_molt,nymphs,nymph_molt,adults)
  
  for (each_day in days){
    
    #Get current week
    week_in <- week(as.Date(each_day,structure(0, class = "Date")))
    
    #get current week temp
    exp_temp <- temp_data$stand_temp[week_in]
    exp_prec <- precip_data$stand_prec[week_in]
    
    #Calculated environmental impacts on molt
    
    p_larv_molt <- br_larvae_molt + (exp_temp_param_larv_molt * exp_temp) + (exp_precip_param_larv_molt * exp_prec)
    
    #Don't let this go below 0! Doesn't make sense biologically!
    if (p_larv_molt < 0){
      p_larv_molt = 0
    } else if (p_larv_molt > max_larvae_molt){
      p_larv_molt = max_larvae_molt
    }
    
    p_nymph_molt <- br_nymph_molt + (exp_temp_param_nymph_molt * exp_temp) + (exp_precip_param_nymph_molt * exp_prec)
    #Don't let this go below 0! Doesn't make sense biologically!
    if (p_nymph_molt < 0){
      p_nymph_molt = 0
    } else if (p_nymph_molt > max_nymph_molt){
      p_nymph_molt = max_nymph_molt
    }
    
    p_adult_breed <- br_adult_breed + (exp_temp_param_adult_breed * exp_temp)
    #Don't let this go below 0! Doesn't make sense biologically!
    if (p_adult_breed < 0){
      p_adult_breed = 0
    } else if (p_adult_breed > .1){
      p_nymph_molt = .1
    }
    
    current_eggs<- tick_record[1,each_day]
    current_larvae <- tick_record[2,each_day]
    current_larvae_molting <- tick_record[3,each_day]
    current_nymphs <- tick_record[4,each_day]
    current_nymphs_molting <- tick_record[5,each_day]
    current_adults <- tick_record[6,each_day]
    
    #We will do some calculations at the start to get rid of redundant terms
    hatching_eggs <- current_eggs*.02
    larvae_to_molt <- current_larvae*p_larv_molt#(4% of larvae will molt)
    new_nymphs <- current_larvae_molting * .04
    nymphs_to_molt <- current_nymphs*p_nymph_molt #(4% of nymphs will molt)
    new_adults <- current_nymphs_molting * .025
    adults_to_breed <- current_adults*p_adult_breed #(30% of adults will breed)
    
    #Here, larvae at t+1 = Survival rate (20 day survival) - outflow to molt + influx from adults laying eggs
    larvae_updated <- current_larvae*.95 - larvae_to_molt + hatching_eggs
    
    #Larvae molting at t+1 = Survival rate (Very high) + influx from larvae finding hosts - new nymphs
    ln_molting_updated <- current_larvae_molting + larvae_to_molt - new_nymphs
    
    #Nymphs at t+1 = Survival rate (40 day survival) + influx of new molts - nymphs that found hosts
    nymphs_updated <- current_nymphs*.975 + new_nymphs - nymphs_to_molt
    
    #Nymphs molting t+1 = Survival rate (Very high) + influx of nymphs - molts to adults
    na_molting_updated <- current_nymphs_molting  + nymphs_to_molt - new_adults
    
    #Adults = survival rate of adults (200 day survival) + influx of nymphs - outflow of breeding adults
    adults_updated <- current_adults*.995 + new_adults - adults_to_breed
    
    #Adults -> eggs 
    eggs_updated <- current_eggs + 2500*adults_to_breed - hatching_eggs
    
    tick_record[1:6,(each_day+1)] <- c(eggs_updated,larvae_updated,ln_molting_updated,nymphs_updated,
                                       na_molting_updated,adults_updated)
  }
  return(tick_record)
}


# Set up the parameters - you can change these!
br_larvae_molt <- .000
br_nymph_molt <- .000

max_larvae_molt <- .01
max_nymph_molt <- .01

exp_temp_param_larv_molt <- 0.0001
exp_temp_param_nymph_molt <- 0.0001
exp_precip_param_larv_molt <- 0.00001
exp_precip_param_nymph_molt <- 0.00001

#Run the simulation
sim_run <- turkey_tick_model_v1(temp_data_avg,precip_data_avg,
                                br_larvae_molt,br_nymph_molt,
                                max_larvae_molt,max_nymph_molt,
                                exp_temp_param_larv_molt,exp_temp_param_nymph_molt,
                                exp_precip_param_larv_molt,exp_precip_param_nymph_molt)

#Plot our data, versus target data
plot(sim_run[2,]/10,type="l",col="forestgreen",lwd=3,lty=1,ylim=c(0,200))
points(sim_run[4,],type="l",col="dodgerblue4",lwd=3,lty=2)
points(sim_run[6,]*10,type="l",col="orchid",lwd=3,lty=3)
points(temp_data_avg$week*7-3.5,temp_data_avg$avg_week_temp,type = "l")
points(ticks_target$mmwrWeek*7-3.5,ticks_target$amblyomma_americanum,cex=1,pch=19)


