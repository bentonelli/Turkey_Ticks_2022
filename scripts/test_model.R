#Test model script

# Setting up the # days the simulation runs for
days <- 1:200

#Read in our data
#Arrays (200 values, represent precipitation, temperature, host phenology, etc)

#Starting population sizes
larvae <- 5000 #by life stage
nymphs <- 50
adults <- 50

#Setting up record keeping
tick_record <- matrix(NA, nrow=3,ncol=length(days)+1)
tick_record[1:3,1] <- c(larvae,nymphs,adults)

for (each_day in days){
  
  print(each_day)
  
  current_larvae <- tick_record[1,each_day]
  current_nymphs <- tick_record[2,each_day]
  current_adults <- tick_record[3,each_day]
  
  #We will say that X% of larvae each day die :(

  #We will say 1% find a host and become nymphs
  # Example with other data larvae_death_rate <- prec[each_day]*.1 + temp[each_day]*.2
  
  #So .98 indicates loss to nymph compartment
  larvae_updated <- current_larvae*.98 - current_larvae*.02
  
  #New nymphs molting
  
  #We will say X% of our nymphs die each day again, and then some become adults (3%)
  nymphs_updated <- current_nymphs*.99 + current_larvae*.02 - current_nymphs*.01
  
  #New adults molting
  
  #Finally, update adult compartment
  adults_updated <- current_adults*.995 + current_nymphs*.01
  
  #Adults -> Larvae (new compartment)
  
  tick_record[1:3,(each_day+1)] <- c(larvae_updated,nymphs_updated,adults_updated)
}

plot(tick_record[1,],type="l",col="forestgreen",lwd=3,lty=1)
points(tick_record[2,],type="l",col="dodgerblue4",lwd=3,lty=2)
points(tick_record[3,],type="l",col="orchid",lwd=3,lty=3)
