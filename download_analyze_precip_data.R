# Code to analyze precipitation data near phenocam sites

# Set working directory
setwd("~/hls-phenocam-plots")

# Libraries
library(dplyr)
library(ggplot2)
library(lubridate)

# Import data from CSV
#
# standard_precip_through_2022.csv is monthly precip downloaded from EDI
monthly_JER_precip <- read.csv('data/precip/standard_precip_through_2022.csv', header=TRUE)
monthly_JER_precip <- select(monthly_JER_precip,c(name,pasture,year,month,sample_date,prec_in,notes))
monthly_JER_precip$sample_date <- as.Date(monthly_JER_precip$sample_date, '%m/%d/%Y')
monthly_JER_precip$prec_in <- monthly_JER_precip$prec_in*25.4 #convert inches to mm
monthly_JER_precip <- rename(monthly_JER_precip, prec_mm = prec_in)
#
# lter_standard_precip_through_2022.csv is precip collected at LTER gauges in
# the Jornada at approximately monthly intervals
standard_LTER_precip <- read.csv('data/precip/lter_standard_precip_through_2022.csv', header=TRUE)
standard_LTER_precip <- select(standard_LTER_precip,c(date,site,ppt_mm,qflag,comment))
standard_LTER_precip$date <- as.Date(standard_LTER_precip$date, '%Y-%m-%d')

# Aggregate LTER standard precip into monthly precip
monthly_LTER_precip <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("name", "year", "month", "prec_mm")
colnames(monthly_LTER_precip) <- x
rm(x)
for (gauge in unique(standard_LTER_precip$site)){
  for (year in 1989:2022){
    for (month in 1:12){
      subset <- standard_LTER_precip[standard_LTER_precip$site==gauge,]
      subset <- subset[order(subset$date),]
      # Find nearest start date for the month
      start_of_current_month <- make_date(year=year, month=month, day=1)
      start_of_current_month <- subset[which.min(abs(start_of_current_month-subset$date)),]$date
      # Find nearest end date for the month
      if (month == 12){
        end_of_current_month <- make_date(year=year+1, month=1, day=1)
      }
      else {
        end_of_current_month <- make_date(year=year, month=month+1, day=1)
      }
      end_of_current_month <- subset[which.min(abs(end_of_current_month-subset$date)),]$date
      # Calculate approximate monthly precip
      monthly_precip_mm <- sum(subset[subset$date > start_of_current_month &
                                      subset$date <= end_of_current_month,]$ppt_mm, na.rm=TRUE)
      new_row <- data.frame(name=gauge, year=year, month=month, prec_mm=monthly_precip_mm)
      monthly_LTER_precip <- rbind(monthly_LTER_precip,new_row)
    }
  }
}
rm(gauge,year,month,subset,start_of_current_month,end_of_current_month,monthly_precip_mm,new_row)

# Standard precip aggregation
monthly_precip <- rbind(select(monthly_JER_precip,c(name,year,month,prec_mm)),monthly_LTER_precip)
monthly_precip <- monthly_precip[monthly_precip$year >= 1989,]
annual_site_precip <- data.frame(matrix(ncol = 10, nrow = 0))
x <- c("site", "year", "total_rainfall_annual", "total_rainfall_winter", "total_rainfall_winter_spring", "total_rainfall_spring", "total_rainfall_spring_summer", "total_rainfall_monsoon", "total_rainfall_summer_fall", "total_rainfall_fall")
colnames(annual_site_precip) <- x
rm(x)

for (gauge in unique(monthly_precip$name)){
  for (year in unique(monthly_precip[monthly_precip$name==gauge,]$year)){
    subset <- monthly_precip[monthly_precip$name == gauge & monthly_precip$year == year,]
    # Calculate total precipitation by calendar year
    if (nrow(subset[subset$year == year,]) != 12){
      total_precip_calendar <- NA
    }
    else {
      total_precip_calendar <- sum(subset[subset$year == year,]$prec_mm, na.rm = TRUE)
    }
    # Calculate total winter precipitation (January-March)
    total_precip_winter <- subset[subset$year == year,]
    total_precip_winter <- total_precip_winter[total_precip_winter$month <= 3,]
    if (nrow(total_precip_winter) != 3){
      total_precip_winter <- NA
    }
    else {
      total_precip_winter <- sum(total_precip_winter$prec_mm, na.rm = TRUE)
    }
    # Calculate total spring precipitation (April-June)
    total_precip_spring <- subset[subset$year == year,]
    total_precip_spring <- total_precip_spring[total_precip_spring$month >= 4 & total_precip_spring$month <= 6,]
    if (nrow(total_precip_spring) != 3){
      total_precip_spring <- NA
    }
    else {
      total_precip_spring <- sum(total_precip_spring$prec_mm, na.rm = TRUE)
    }
    # Calculate total precipitation during monsoon season (July, August, September)
    total_precip_monsoon <- subset[subset$year == year,]
    total_precip_monsoon <- total_precip_monsoon[total_precip_monsoon$month >= 7 & total_precip_monsoon$month <= 9,]
    if (nrow(total_precip_monsoon) != 3){
      total_precip_monsoon <- NA
    }
    else {
      total_precip_monsoon <- sum(total_precip_monsoon$prec_mm, na.rm = TRUE)
    }
    # Calculate total fall precipitation (October-December)
    total_precip_fall <- subset[subset$year == year,]
    total_precip_fall <- total_precip_fall[total_precip_fall$month >= 10,]
    if (nrow(total_precip_fall) != 3){
      total_precip_fall <- NA
    }
    else {
      total_precip_fall <- sum(total_precip_fall$prec_mm, na.rm = TRUE)
    }
    
    # Populate dataframe
    annual_site_precip <- rbind(annual_site_precip,
                                data.frame(site=gauge,
                                      year=year,
                                      total_rainfall_annual=total_precip_calendar,
                                      total_rainfall_winter=total_precip_winter,
                                      total_rainfall_winter_spring=total_precip_winter+total_precip_spring,
                                      total_rainfall_spring=total_precip_spring,
                                      total_rainfall_spring_summer=total_precip_spring+total_precip_monsoon,
                                      total_rainfall_monsoon=total_precip_monsoon,
                                      total_rainfall_summer_fall=total_precip_monsoon+total_precip_fall,
                                      total_rainfall_fall=total_precip_fall))
  }
}
rm(subset,gauge,total_precip_calendar,total_precip_monsoon,total_precip_winter,total_precip_spring,year)

# Annual site precip for selected gauges and years
annual_site_precip_selected_gauges <- annual_site_precip[annual_site_precip$site == "ARISTIDA" |
                                                           #annual_site_precip$site == "DONA ANA" |
                                                           #annual_site_precip$site == "IBPE" |
                                                           annual_site_precip$site == "SAND" |
                                                           annual_site_precip$site == "RAGGED TANK" |
                                                           annual_site_precip$site == "RABBIT" |
                                                           annual_site_precip$site == "HEADQUARTERS" |
                                                           annual_site_precip$site == "WEST WELL" |
                                                           annual_site_precip$site == "CO-OP WELL",]

annual_site_precip_selected_years <- annual_site_precip_selected_gauges[annual_site_precip_selected_gauges$year > 2013,]

# Plot precipitation totals over different seasons

# Plot annual precip for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=total_rainfall_annual, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("Annual precipitation (mm)") +
  ggtitle("Annual precipitation")
p

# Plot winter precip for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=total_rainfall_winter, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("Winter precipitation (mm)") +
  ggtitle("Winter precipitation (Jan-Mar)")
p

# Plot winter and spring precip for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=total_rainfall_winter_spring, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("Winter and spring precipitation (mm)") +
  ggtitle("Winter and spring precipitation (Jan-Jun)")
p

# Plot spring precip for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=total_rainfall_spring, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("Spring precipitation (mm)") +
  ggtitle("Spring precipitation (Apr-Jun)")
p

# Plot spring and summer precip for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=total_rainfall_spring_summer, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("Spring and summer precipitation (mm)") +
  ggtitle("Spring and summer precipitation (Apr-Sep)")
p

# Plot monsoon season precip for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=total_rainfall_monsoon, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("Monsoon season precipitation (mm)") +
  ggtitle("Monsoon season precipitation (Jul-Sep)")
p

# Plot summer and fall precip for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=total_rainfall_summer_fall, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("Summer and fall precipitation (mm)") +
  ggtitle("Summer and fall precipitation (Jul-Dec)")
p

# Calculate SMD

annual_site_precip_selected_years <- cbind(annual_site_precip_selected_years,
                                           smd_annual = 0,
                                           smd_winter = 0,
                                           smd_winter_spring = 0,
                                           smd_spring = 0,
                                           smd_spring_summer = 0,
                                           smd_monsoon = 0,
                                           smd_summer_fall = 0,
                                           smd_fall = 0)

# Calculate long-term (1989-2022) SMD for each gauge
for(gauge in unique(annual_site_precip_selected_gauges$site)){
  subset <- annual_site_precip_selected_gauges[annual_site_precip_selected_gauges$site==gauge,]
  # Calendar year
  long_term_mean_precip <- mean(subset[subset$site==gauge,]$total_rainfall_annual)
  long_term_sdev_precip <- sd(subset[subset$site==gauge,]$total_rainfall_annual)
  # Winter
  long_term_mean_winter_precip <- mean(subset[subset$site==gauge,]$total_rainfall_winter)
  long_term_sdev_winter_precip <- sd(subset[subset$site==gauge,]$total_rainfall_winter)
  # Winter and spring
  long_term_mean_winter_spring_precip <- mean(subset[subset$site==gauge,]$total_rainfall_winter_spring)
  long_term_sdev_winter_spring_precip <- sd(subset[subset$site==gauge,]$total_rainfall_winter_spring)
  # Spring
  long_term_mean_spring_precip <- mean(subset[subset$site==gauge,]$total_rainfall_spring)
  long_term_sdev_spring_precip <- sd(subset[subset$site==gauge,]$total_rainfall_spring)
  # Spring and summer
  long_term_mean_spring_summer_precip <- mean(subset[subset$site==gauge,]$total_rainfall_spring_summer)
  long_term_sdev_spring_summer_precip <- sd(subset[subset$site==gauge,]$total_rainfall_spring_summer)
  # Monsoon (summer)
  long_term_mean_monsoon_precip <- mean(subset[subset$site==gauge,]$total_rainfall_monsoon)
  long_term_sdev_monsoon_precip <- sd(subset[subset$site==gauge,]$total_rainfall_monsoon)
  # Summer and fall
  long_term_mean_summer_fall_precip <- mean(subset[subset$site==gauge,]$total_rainfall_summer_fall)
  long_term_sdev_summer_fall_precip <- sd(subset[subset$site==gauge,]$total_rainfall_summer_fall)
  # Fall
  long_term_mean_fall_precip <- mean(subset[subset$site==gauge,]$total_rainfall_fall)
  long_term_sdev_fall_precip <- sd(subset[subset$site==gauge,]$total_rainfall_fall)
  
  subset <- subset[subset$year > 2013,]
  subset$smd_annual <- (subset$total_rainfall_annual - long_term_mean_precip) / long_term_sdev_precip
  subset$smd_winter <- (subset$total_rainfall_winter - long_term_mean_winter_precip) / long_term_sdev_winter_precip
  subset$smd_winter_spring <- (subset$total_rainfall_winter_spring - long_term_mean_winter_spring_precip) / long_term_sdev_winter_spring_precip
  subset$smd_spring <- (subset$total_rainfall_spring - long_term_mean_spring_precip) / long_term_sdev_spring_precip
  subset$smd_spring_summer <- (subset$total_rainfall_spring_summer - long_term_mean_spring_summer_precip) / long_term_sdev_spring_summer_precip
  subset$smd_monsoon <- (subset$total_rainfall_monsoon - long_term_mean_monsoon_precip) / long_term_sdev_monsoon_precip
  subset$smd_summer_fall <- (subset$total_rainfall_summer_fall - long_term_mean_summer_fall_precip) / long_term_sdev_summer_fall_precip
  subset$smd_fall <- (subset$total_rainfall_fall - long_term_mean_fall_precip) / long_term_sdev_fall_precip
  
  annual_site_precip_selected_years[annual_site_precip_selected_years$site==gauge,]$smd_annual <- subset$smd_annual
  annual_site_precip_selected_years[annual_site_precip_selected_years$site==gauge,]$smd_winter <- subset$smd_winter
  annual_site_precip_selected_years[annual_site_precip_selected_years$site==gauge,]$smd_winter_spring <- subset$smd_winter_spring
  annual_site_precip_selected_years[annual_site_precip_selected_years$site==gauge,]$smd_spring <- subset$smd_spring
  annual_site_precip_selected_years[annual_site_precip_selected_years$site==gauge,]$smd_spring_summer <- subset$smd_spring_summer
  annual_site_precip_selected_years[annual_site_precip_selected_years$site==gauge,]$smd_monsoon <- subset$smd_monsoon
  annual_site_precip_selected_years[annual_site_precip_selected_years$site==gauge,]$smd_summer_fall <- subset$smd_summer_fall
  annual_site_precip_selected_years[annual_site_precip_selected_years$site==gauge,]$smd_fall <- subset$smd_fall
}

# Calculate mean and sdev of spring-summer rainfall
long_term_mean_spring_summer_all_gauges <- mean(annual_site_precip_selected_gauges$total_rainfall_spring_summer)
long_term_sdev_spring_summer_all_gauges <- sd(annual_site_precip_selected_gauges$total_rainfall_spring_summer)
upper_thresh_spring_summer <- long_term_mean_spring_summer_all_gauges + 0.5*long_term_sdev_spring_summer_all_gauges
lower_thresh_spring_summer <- long_term_mean_spring_summer_all_gauges - 0.5*long_term_sdev_spring_summer_all_gauges
annual_site_precip_selected_years <- cbind(annual_site_precip_selected_years, wet_dry_ave_spring_summer = "average")
annual_site_precip_selected_years[annual_site_precip_selected_years$total_rainfall_spring_summer>=upper_thresh_spring_summer,]$wet_dry_ave_spring_summer <- "wet"
annual_site_precip_selected_years[annual_site_precip_selected_years$total_rainfall_spring_summer<=lower_thresh_spring_summer,]$wet_dry_ave_spring_summer <- "dry"

# Calculate mean and sdev of winter rainfall
long_term_mean_winter_all_gauges <- mean(annual_site_precip_selected_gauges$total_rainfall_winter)
long_term_sdev_winter_all_gauges <- sd(annual_site_precip_selected_gauges$total_rainfall_winter)
upper_thresh_winter <- long_term_mean_winter_all_gauges + 0.5*long_term_sdev_winter_all_gauges
lower_thresh_winter <- long_term_mean_winter_all_gauges - 0.5*long_term_sdev_winter_all_gauges
annual_site_precip_selected_years <- cbind(annual_site_precip_selected_years, wet_dry_ave_winter = "average")
annual_site_precip_selected_years[annual_site_precip_selected_years$total_rainfall_winter>=upper_thresh_winter,]$wet_dry_ave_winter <- "wet"
annual_site_precip_selected_years[annual_site_precip_selected_years$total_rainfall_winter<=lower_thresh_winter,]$wet_dry_ave_winter <- "dry"

# Plot SMD over different seasons

# Plot annual SMD for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=smd_annual, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("SMD (1989-2022 annual precipitation)") +
  ggtitle("Annual precipitation SMD")
p

# Plot winter SMD for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=smd_winter, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("SMD (1989-2022 winter precipitation)") +
  ggtitle("Winter precipitation SMD (Jan-Mar)")
p

# Plot winter and spring SMD for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=smd_winter_spring, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("SMD (1989-2022 winter precipitation)") +
  ggtitle("Winter and spring precipitation SMD (Jan-Jun)")
p

# Plot spring SMD for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=smd_spring, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("SMD (1989-2022 spring precipitation)") +
  ggtitle("Spring precipitation SMD (Apr-Jun)")
p

# Plot spring and summer SMD for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=smd_spring_summer, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("SMD (1989-2022 spring and summer precipitation)") +
  ggtitle("Spring and summer precipitation SMD (Apr-Sep)")
p

# Plot monsoon season SMD for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=smd_monsoon, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("SMD (1989-2022 monsoon season precipitation)") +
  ggtitle("Monsoon season precipitation SMD (Jul-Sep)")
p

# Plot summer and fall SMD for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=smd_summer_fall, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("SMD (1989-2022 summer and fall precipitation)") +
  ggtitle("Summer and fall precipitation SMD (Jul-Dec)")
p

# Plot fall SMD for each gauge
p <- ggplot(data = annual_site_precip_selected_years, aes(x=year, y=smd_fall, color=site)) +
  geom_point() +
  geom_line() +
  xlab("Calendar year") +
  ylab("SMD (1989-2022 fall precipitation)") +
  ggtitle("Fall precipitation SMD (Oct-Dec)")
p

# Use SMD to determine wet, dry, and average years

wet_dry_average <- function(value, threshold){
  if(value >= threshold){
    rainfall = "wet"
  }
  else if(value <= -1*threshold){
    rainfall = "dry"
  }
  else{
    rainfall = "average"
  }
  rainfall
}

# 0.5 threshold
thresh = 0.5
annual_site_precip_selected_years <- cbind(annual_site_precip_selected_years,
                                           rainfall_annual = sapply(annual_site_precip_selected_years$smd_annual,wet_dry_average,thresh),
                                           rainfall_winter = sapply(annual_site_precip_selected_years$smd_winter,wet_dry_average,thresh),
                                           rainfall_winter_spring = sapply(annual_site_precip_selected_years$smd_winter_spring,wet_dry_average,thresh),
                                           rainfall_spring = sapply(annual_site_precip_selected_years$smd_spring,wet_dry_average,thresh),
                                           rainfall_spring_summer = sapply(annual_site_precip_selected_years$smd_spring_summer,wet_dry_average,thresh),
                                           rainfall_monsoon = sapply(annual_site_precip_selected_years$smd_monsoon,wet_dry_average,thresh),
                                           rainfall_summer_fall = sapply(annual_site_precip_selected_years$smd_summer_fall,wet_dry_average,thresh),
                                           rainfall_fall = sapply(annual_site_precip_selected_years$smd_fall,wet_dry_average,thresh)
                                           )

# Make a dataframe with phenocam site and corresponding rain gauge
phenocam_rain_gauge_pairs <- rbind(data.frame(phenocam_name = "jershrubland", site = "ARISTIDA"),
                                   data.frame(phenocam_name = "jershrubland2", site = "ARISTIDA"),
                                   data.frame(phenocam_name = "jernovel", site = "ARISTIDA"),
                                   data.frame(phenocam_name = "jernovel2", site = "ARISTIDA"),
                                   data.frame(phenocam_name = "jergrassland", site = "WEST WELL"),
                                   data.frame(phenocam_name = "jergrassland2", site = "CO-OP WELL"),
                                   data.frame(phenocam_name = "jerbajada", site = "RAGGED TANK"),
                                   data.frame(phenocam_name = "jernort", site = "RABBIT"),
                                   data.frame(phenocam_name = "ibp", site = "WEST WELL"),
                                   data.frame(phenocam_name = "jernwern", site = "HEADQUARTERS"),
                                   data.frame(phenocam_name = "NEON.D14.JORN.DP1.00033", site = "WEST WELL"),
                                   data.frame(phenocam_name = "jersand", site = "SAND"))

# Save results
save(list=c("annual_site_precip_selected_years","phenocam_rain_gauge_pairs"), file="outputs/annual_precip.RData")
