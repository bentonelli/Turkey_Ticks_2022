### Script to process eBird data ####

library(auk)
library(lubridate)
library(dplyr)
filenames = list.files("data/ebird_count_tick_forecasting/txt_files/",pattern = "*.txt")
file_rec <- c()
for (each_ebird_file in filenames){
  print(each_ebird_file)
  #read in checklist, eliminating duplicate group checlists and rolling up taxonomy
  f_ebd <- read_ebd(paste("data/ebird_count_tick_forecasting/txt_files/",each_ebird_file,sep=""))
  zf_checklists <- as.data.frame(table(f_ebd$checklist_id,f_ebd$common_name))
  colnames(zf_checklists) <- c("checklist_id","common_name","p_a")
  #Add date back in
  checklist_dates <- select(f_ebd,c("checklist_id","observation_date"))
  checklist_dates <- checklist_dates[!duplicated(checklist_dates),]
  
  zf_checklists <- merge(zf_checklists,checklist_dates)
  turkey_checklists <- zf_checklists %>%
    group_by(week(observation_date)) %>%
    filter(common_name == "Wild Turkey") %>%
    summarise(perc_checklists=sum(p_a)/n())
  
  colnames(turkey_checklists) <- c("week","perc_checklists")
  turkey_checklists$county <- f_ebd$county[1]
  file_rec <- rbind(file_rec,turkey_checklists)
}

write.csv(file_rec,"turkey_trends_week.csv")


