#Turkey tick mass simulator
# Function to run tick population model LOTS OF TIMES, V1.0 
library(dplyr)
library(lubridate)
library(truncnorm) #New package for most people

### DATA ####

#Read in our training data - this is site specific, feel free to change the site
ticks_target <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)
ticks_target <- filter(ticks_target,siteID == "SCBI")

#Read in our temperature data - again, site specific
temp_data <- read.csv("https://raw.githubusercontent.com/bentonelli/Turkey_Ticks_2022/master/data/temp_data/Weekly%20Average%20Temperature/%20SCBI%20.csv")
temp_data$year <- as.numeric(substr(temp_data$Week, 1, 4)) 
temp_data$week_num <- as.numeric(substr(temp_data$Week, 6, 7)) 

#Process so that it fits a workable format
temp_data_avg <- temp_data %>% 
  group_by(week=temp_data$week_num) %>%
  summarise(avg_week_temp = mean(tempMean))
temp_data_avg$stand_temp <- (temp_data_avg$avg_week_temp - mean(temp_data_avg$avg_week_temp))


#Read in our temperature data - again, site specific
precip_data <- read.csv("https://raw.githubusercontent.com/bentonelli/Turkey_Ticks_2022/master/data/Combined_weekly_precip_data.csv")
precip_data <- filter(precip_data,siteID == "SCBI")

#Process so that it fits a workable format
precip_data_avg <- precip_data %>%
  group_by(week) %>%
  summarise(mean_week_precip = mean(sum_weekly_total_precip_bySite))
precip_data_avg$stand_prec <- precip_data_avg$mean_week_precip - mean(precip_data_avg$mean_week_precip)

##deerdata####
birth_rate <- 0.66
initial_pop = 100
n_fawns <- birth_rate*initial_pop
n_fawns
m <- 19
stdv <-2
d_surv <- 0.99 
birth_dates <- round(rnorm(n_fawns,m,stdv))
birth_dates
pop_rec <- c(initial_pop)
for(each_week in 1:52){
  num_births = sum(birth_dates==each_week)
  d_pop_1 = pop_rec[each_week]*d_surv + num_births 
  pop_rec <- c(pop_rec, d_pop_1)
}
pop_rec - mean(pop_rec)
pop_rec_stand <- pop_rec - mean(pop_rec)


### Model Code ####

turkey_tick_model_v1 <- function(temp_data = temp_data_avg,precip_data = precip_data_avg,deer_data=pop_rec_stand,
                                 br_larvae_molt=.001,br_nymph_molt=.001,br_adult_breed=.001,
                                 max_larvae_molt=.01,max_nymph_molt=.01,max_adult_breed=.01,
                                 exp_temp_param_larv_molt=0.0000,exp_temp_param_nymph_molt = 0.0001,
                                 exp_precip_param_larv_molt=0.00001,exp_precip_param_nymph_molt = 0.00001,
                                 exp_deer_param_larv_molt=0.00001,
                                 exp_deer_param_nymph_molt = 0.0001,
                                 exp_deer_param_adult_breed= 0.00001,
                                 exp_temp_param_adult_breed=0.0002,
                                 start_date=1,starting_nymphs = 10,starting_adults=10){
  # Setting up the # days the simulation runs for
  days <- start_date:365
  #Starting population sizes
  eggs <- 0
  larvae <- 0 #by life stage
  larvae_molt <- 0
  nymphs <- starting_nymphs
  nymph_molt <- 0
  adults <- starting_adults
  
  #Setting up record keeping
  tick_record <- matrix(NA, nrow=6,ncol=366)
  tick_record[1:6,start_date] <- c(eggs,larvae,larvae_molt,nymphs,nymph_molt,adults)
  
  for (each_day in days){
    
    #Get current week
    week_in <- week(as.Date(each_day,structure(0, class = "Date")))
    
    #get current week temp
    exp_temp <- temp_data$stand_temp[week_in]
    exp_prec <- precip_data$stand_prec[week_in]
    exp_deer <- deer_data[week_in]
    #Calculated environmental impacts on molt
    
    p_larv_molt <- br_larvae_molt + (exp_temp_param_larv_molt * exp_temp) + (exp_precip_param_larv_molt * exp_prec) + (exp_deer_param_larv_molt*exp_deer)
    #Don't let this go below 0! Doesn't make sense biologically!
    if (p_larv_molt < 0){
      p_larv_molt = 0
    } else if (p_larv_molt > max_larvae_molt){
      p_larv_molt = max_larvae_molt
    }
    
    p_nymph_molt <- br_nymph_molt + (exp_temp_param_nymph_molt * exp_temp) + (exp_precip_param_nymph_molt * exp_prec) + (exp_deer_param_nymph_molt*exp_deer)
    #Don't let this go below 0! Doesn't make sense biologically!
    if (p_nymph_molt < 0){
      p_nymph_molt = 0
    } else if (p_nymph_molt > max_nymph_molt){
      p_nymph_molt = max_nymph_molt
    }
    
    p_adult_breed <- br_adult_breed + (exp_temp_param_adult_breed * exp_temp) + (exp_deer_param_adult_breed*exp_deer)
    #Don't let this go below 0! Doesn't make sense biologically!
    if (p_adult_breed < 0){
      p_adult_breed = 0
    } else if (p_adult_breed > max_adult_breed){
      p_adult_breed = max_adult_breed
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
    nymphs_updated <- current_nymphs*.95 + new_nymphs - nymphs_to_molt #.975
    
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

# ### Single Model Run ####
# 
# br_larvae_molt <- .002
# br_nymph_molt <- .001
# br_adult_breed <- .001
# 
# max_larvae_molt <- .2
# max_nymph_molt <- .2
# max_adult_breed <-.2
# 
# exp_temp_param_larv_molt <- 0.002
# exp_temp_param_nymph_molt <- 0.0001
# exp_precip_param_larv_molt <- 0.00001
# exp_precip_param_nymph_molt <- 0.00001
# exp_temp_param_adult_breed <- .00005
# 
# start_date <- 1
# starting_nymphs <- 10
# starting_adults <- 20
# 
# sim_run <- turkey_tick_model_v1(temp_data = temp_data_avg,precip_data = precip_data_avg,
#                                 br_larvae_molt = br_larvae_molt,br_nymph_molt = br_nymph_molt,br_adult_breed = br_adult_breed,
#                                 max_larvae_molt = max_larvae_molt,max_nymph_molt = max_nymph_molt,max_adult_breed = max_adult_breed,
#                                 exp_temp_param_larv_molt = exp_temp_param_larv_molt,exp_temp_param_nymph_molt = exp_temp_param_nymph_molt,
#                                 exp_precip_param_larv_molt = exp_precip_param_larv_molt,exp_precip_param_nymph_molt = exp_precip_param_nymph_molt,
#                                 exp_temp_param_adult_breed = exp_temp_param_adult_breed,
#                                 start_date=start_date,starting_nymphs = starting_nymphs,starting_adults=starting_adults)
# 
# plot(sim_run[2,]/10,type="l",col="forestgreen",lwd=3,lty=1,ylim=c(0,200))
# points(sim_run[4,],type="l",col="dodgerblue4",lwd=3,lty=2)
# points(sim_run[6,]*10,type="l",col="orchid",lwd=3,lty=3)
# points(temp_data_avg$week*7-3.5,temp_data_avg$avg_week_temp,type = "l")
# points(60+ticks_target$mmwrWeek*7-3.5,ticks_target$amblyomma_americanum,cex=1,pch=19)




### Multi Model Run ####

# Get target data
week_df <- as.data.frame(1:52)
colnames(week_df) <- "week"
target_data <- ticks_target %>% 
  group_by(week = mmwrWeek) %>%
  summarise(mean_ticks = mean(amblyomma_americanum))
target_data <- merge(week_df,target_data,all = TRUE)
target_data$mean_ticks[which(is.na(target_data$mean_ticks))] <- 0

#Run a bunch of simulations
num_sims <- 100
best_err_cutoff <- 9999999999
param_set_rec <- c()
summa <- vector(mode = "list", length = 20)

for (m in 1:20){
for (n in 1:num_sims){
  # print(n)
  
  #Draw random parameters - you need to set reasonable ranges here.
  
  br_larvae_molt <- rnorm(1,0.003,.0001)
  br_nymph_molt <- rnorm(1,0.01,.001)
  br_adult_breed <- rnorm(1,0,.001)
  
  max_larvae_molt <- .08
  max_nymph_molt <- .1
  max_adult_breed <- .4
  
  exp_temp_param_larv_molt <- rtruncnorm(1,a=0,mean=.001,sd=.0001)
  exp_temp_param_nymph_molt <- rtruncnorm(1,a=0.25,mean=.5,sd=.01)
  exp_precip_param_larv_molt <- rnorm(1,mean=.0001,sd=.00003)
  exp_precip_param_nymph_molt <- rnorm(1,mean=.002,sd=.0005)
  exp_temp_param_adult_breed <- rnorm(1,mean=.0005,sd=.00005)
  exp_deer_param_nymph_molt <- rnorm(1,mean=.08,sd=.001)
  exp_deer_param_larv_molt <- rnorm(1,mean=.004,sd=.0005)
  exp_deer_param_adult_breed <- rnorm(1,mean=.0005,sd=.0005)
  
  nymph_start <- rtruncnorm(1,a=0,mean=100,sd=5)
  adult_start <- rtruncnorm(1,a=0,mean=50,sd=5)
  
  start_date <- 20 ## Keep this at 1 for now
  
  #Run simulation
  sim_run <- turkey_tick_model_v1(temp_data = temp_data_avg,precip_data = precip_data_avg,
                                  #Base rate of molting, breeding
                                  br_larvae_molt = br_larvae_molt,
                                  br_nymph_molt = br_nymph_molt,
                                  br_adult_breed = br_adult_breed,
                                  #Max rates of molting, breeding
                                  max_larvae_molt = max_larvae_molt,
                                  max_nymph_molt = max_nymph_molt,
                                  max_adult_breed = max_adult_breed,
                                  #Temperature
                                  exp_temp_param_larv_molt = exp_temp_param_larv_molt,
                                  exp_temp_param_nymph_molt = exp_temp_param_nymph_molt,
                                  exp_temp_param_adult_breed = exp_temp_param_adult_breed, 
                                  #Precip.
                                  exp_precip_param_larv_molt = exp_precip_param_larv_molt,
                                  exp_precip_param_nymph_molt = exp_precip_param_nymph_molt,
                                  #deer 
                                  exp_deer_param_larv_molt = exp_deer_param_larv_molt,
                                  exp_deer_param_nymph_molt = exp_deer_param_nymph_molt,
                                  exp_deer_param_adult_breed = exp_deer_param_adult_breed,
                                  #Start conditions
                                  start_date=start_date,
                                  starting_nymphs = nymph_start,
                                  starting_adults= adult_start
  )
  
  
  #Compare to real data
  sim_week_data <- sum(abs(sim_run[4,seq(3,365,by=7)] - target_data$mean_ticks),na.rm=TRUE)
  
  # Translation factor - this finds the best fit with the observed data based on 
  # when the simulation starts. We can do this because our current start date of 
  # Jan 1 is ultimately arbitrary
  
  # These factors represent change in WEEKS, such that a value of -1 means we are testing
  # the fit if we started the simulation one week earlier (Dec 24th)
  factor <- seq(-20,20,by=1)
  best_translation <- 99999999999
  trans_week <- -999
  for (each_factor in factor){
    translated_target_data <- rep(NA,52)
    if (each_factor<=0){
      translated_target_data[1:(52+each_factor)] <- abs(sim_run[4,seq(3,365,by=7)])[(1-each_factor):52]
    } else if (each_factor > 0){
      translated_target_data[(1+each_factor):52] <- abs(sim_run[4,seq(3,365,by=7)])[1:(52-each_factor)]
    } else {
      print("error")
    }
    trans_err <- sum(abs(translated_target_data - target_data$mean_ticks),na.rm = TRUE)
    if (trans_err < best_translation){
      trans_week <- each_factor
      best_translation <- trans_err
      trans_tick_pred <- translated_target_data
      total_ticks <- sum(translated_target_data,na.rm=TRUE)
    }
  }
  
  #Save best model for plotting
  if(best_translation < best_err_cutoff){
    best_err <- trans_tick_pred
    best_err_cutoff <- best_translation
  }
  #Save parameters, error
  param_set <- c(best_translation,total_ticks,trans_week,
                 br_larvae_molt,br_nymph_molt,br_adult_breed,
                 max_larvae_molt,max_nymph_molt,max_adult_breed,
                 exp_temp_param_larv_molt,exp_temp_param_nymph_molt,exp_temp_param_adult_breed,
                 exp_precip_param_larv_molt,exp_precip_param_nymph_molt,
                 exp_deer_param_larv_molt,
                 exp_deer_param_nymph_molt,
                 exp_deer_param_adult_breed,
                 start_date,nymph_start,adult_start)
  param_set_rec <- rbind(param_set_rec,param_set)
  
}

### Analysis ####
#See how the best model did!
plot(target_data$week,target_data$mean_ticks,pch=19,col="black")
points(1:52,best_err,type="l",col="dodgerblue4",lwd=3,lty=2)

#Save and view our data
param_set_rec <- as.data.frame(param_set_rec)
colnames(param_set_rec) <- c("best_translation","total_ticks","trans_week",
                             "br_larvae_molt","br_nymph_molt","br_adult_breed",
                             "max_larvae_molt","max_nymph_molt","max_adult_breed",
                             "exp_temp_param_larv_molt","exp_temp_param_nymph_molt","exp_temp_param_adult_breed",
                             "exp_precip_param_larv_molt","exp_precip_param_nymph_molt",
                             "exp_deer_param_larv_molt",
                             "exp_deer_param_nymph_molt",
                             "exp_deer_param_adult_breed",
                             "start_date","nymph_start","adult_start")
# View(param_set_rec)
summa[[m]] <- summary(param_set_rec$best_translation)
}

summa
yey <- param_set_rec[(param_set_rec$best_translation == min(param_set_rec$best_translation)),]

# #Look at relationships between starting parameters and simulation errr
# plot(param_set_rec$nymph_start,param_set_rec$best_translation,cex=.1)
# plot(param_set_rec$adult_start,param_set_rec$best_translation,cex=.1)
# #Molt, breeding base rates
# plot(param_set_rec$br_larvae_molt,param_set_rec$best_translation,cex=.1)
# plot(param_set_rec$br_nymph_molt,param_set_rec$best_translation,cex=.1)
# plot(param_set_rec$br_adult_breed,param_set_rec$best_translation,cex=.1)
# #Temp
# plot(param_set_rec$exp_temp_param_larv_molt,param_set_rec$best_translation,cex=.1)
# plot(param_set_rec$exp_temp_param_nymph_molt,param_set_rec$best_translation,cex=.1)
# plot(param_set_rec$exp_temp_param_adult_breed,param_set_rec$best_translation,cex=.1)
# #Precip
# plot(param_set_rec$exp_precip_param_larv_molt,param_set_rec$best_translation,cex=.1)
# plot(param_set_rec$exp_precip_param_nymph_molt,param_set_rec$best_translation,cex=.1)

