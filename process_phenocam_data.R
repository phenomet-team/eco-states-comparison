## Code to process Phenocam time-series data

## NOTE: Before running this code, please run download_phenocam_data.R or
## load the previously-imported PhenoCam data

# Set working directory
setwd("~/hls-phenocam-plots")

# Load packages
library(dplyr)

# Load previously-imported PhenoCam data
load("outputs/phenocam_imported.RData")
phenocam_data <- select(phenocam_data_3day,date,year,doy,phenocam_name,eco_state,gcc_90,smooth_gcc_90)

# Remove 2023 (not enough data)
phenocam_data <- phenocam_data[phenocam_data$year < 2023,]

############################################################
# Scale the phenocam data (normalize by sensor and by year)
############################################################

# NOTES:
# gcc_scaled_all_years -> GCC scaled over entire time series
# gcc_scaled_yearly -> GCC scaled over each year in the time series

phenocam_data <- cbind(phenocam_data,gcc_scaled_all_years=0,gcc_scaled_yearly=0)

for (camera in unique(phenocam_data$phenocam_name)){
  subset_sensor <- phenocam_data[phenocam_data$phenocam_name==camera,]
  gcc_max <- max(subset_sensor$smooth_gcc_90, na.rm=TRUE)
  gcc_min <- min(subset_sensor$smooth_gcc_90, na.rm=TRUE)
  subset_sensor$gcc_scaled_all_years <- (subset_sensor$smooth_gcc_90 - gcc_min) / (gcc_max - gcc_min)
  phenocam_data[phenocam_data$phenocam_name==camera,]$gcc_scaled_all_years <- subset_sensor$gcc_scaled_all_years
  
  for (year in unique(subset_sensor$year)){
    subset_year <- subset_sensor[subset_sensor$year==year,]
    gcc_max <- max(subset_year$smooth_gcc_90, na.rm=TRUE)
    gcc_min <- min(subset_year$smooth_gcc_90, na.rm=TRUE)
    subset_year$gcc_scaled_yearly <- (subset_year$smooth_gcc_90 - gcc_min) / (gcc_max - gcc_min)
    phenocam_data[phenocam_data$phenocam_name==camera & phenocam_data$year==year,]$gcc_scaled_yearly <- subset_year$gcc_scaled_yearly
    
  }
}

rm(camera,year,subset_sensor,subset_year,gcc_max,gcc_min)

#######################################
# Calculate season start and end dates
#######################################

# NOTES: For thresholding, added the restrictions that SOS must occur after GCC
# reaches its pre-peak minimum value and EOS must occur before GCC reaches its
# post-peak minimum value

# Initialize data frame
season_start_and_end <- data.frame(matrix(ncol = 14, nrow = 0))
x <- c("phenocam_name", "eco_state", "year", "peak", "DOY_min_pre_peak", "DOY_min_post_peak",
       "SOS_10", "SOS_15", "SOS_25", "SOS_50", "EOS_10", "EOS_15", "EOS_25", "EOS_50")
colnames(season_start_and_end) <- x

for (camera in unique(phenocam_data$phenocam_name)){
  for (year in unique(phenocam_data[phenocam_data$phenocam_name==camera,]$year)){
    subset <- phenocam_data[phenocam_data$phenocam_name==camera & phenocam_data$year==year,]
    eco_state <- subset$eco_state[1]
    # Peak - day of the year when GCC is at its maximum value
    peak <- min(subset[subset$gcc_scaled_yearly == 1,]$doy, na.rm = TRUE)
    # Pre-peak minimum - day of the year when GCC is at its minimum value, prior to peak
    min1 <- min(subset[subset$doy < peak,]$gcc_scaled_yearly, na.rm = TRUE)
    min1 <- min(subset[subset$gcc_scaled_yearly == min1,]$doy, na.rm = TRUE)
    # Post-peak minimum - day of the year when GCC is at its minimum value, post peak
    min2 <- min(subset[subset$doy > peak,]$gcc_scaled_yearly, na.rm = TRUE)
    min2 <- max(subset[subset$gcc_scaled_yearly == min2,]$doy, na.rm = TRUE)
    # Start of season - first day of the year when GCC is above the threshold
    SOS_10 <- min(subset[subset$gcc_scaled_yearly > 0.10 & subset$doy >= min1,]$doy, na.rm = TRUE)
    SOS_15 <- min(subset[subset$gcc_scaled_yearly > 0.15 & subset$doy >= min1,]$doy, na.rm = TRUE)
    SOS_25 <- min(subset[subset$gcc_scaled_yearly > 0.25 & subset$doy >= min1,]$doy, na.rm = TRUE)
    SOS_50 <- min(subset[subset$gcc_scaled_yearly > 0.50 & subset$doy >= min1,]$doy, na.rm = TRUE)
    # End of season - last day of the year when GCC is above the threshold
    EOS_10 <- max(subset[subset$gcc_scaled_yearly > 0.10 & subset$doy <= min2,]$doy, na.rm = TRUE)
    EOS_15 <- max(subset[subset$gcc_scaled_yearly > 0.15 & subset$doy <= min2,]$doy, na.rm = TRUE)
    EOS_25 <- max(subset[subset$gcc_scaled_yearly > 0.25 & subset$doy <= min2,]$doy, na.rm = TRUE)
    EOS_50 <- max(subset[subset$gcc_scaled_yearly > 0.50 & subset$doy <= min2,]$doy, na.rm = TRUE)
    new_row <- data.frame(phenocam_name = camera,
                          eco_state = eco_state,
                          year = year,
                          peak = peak,
                          DOY_min_pre_peak = min1,
                          DOY_min_post_peak = min2,
                          SOS_10 = SOS_10,
                          SOS_15 = SOS_15,
                          SOS_25 = SOS_25,
                          SOS_50 = SOS_50,
                          EOS_10 = EOS_10,
                          EOS_15 = EOS_15,
                          EOS_25 = EOS_25,
                          EOS_50 = EOS_50)
    season_start_and_end <- rbind(season_start_and_end,new_row)
  }
}

rm(x,camera,year,subset,eco_state,peak,min1,min2,
   SOS_10,SOS_15,SOS_25,SOS_50,EOS_10,EOS_15,EOS_25,EOS_50,new_row)

# Remove 2023 (not enough data)
season_start_and_end <- season_start_and_end[season_start_and_end$year < 2023,]

###########################
# Calculate cumulative GCC
###########################

# Define function
cumulative_GCC <- function(gcc_df,time_df,start_doy,end_doy,baseline=0){
  df <- data.frame(GCC = gcc_df, DOY = time_df)
  df <- cbind(df, CUMULATIVE_GCC = 0)
  df <- df[df$DOY >= start_doy & df$DOY <= end_doy,]
  df$GCC <- df$GCC - baseline
  # commented-out line below does not work as intended, but is supposed to make
  # it so that if GCC dips below baseline, cumulative GCC does not decrease
  # df[df$GCC < 0,]$GCC <- 0
  df$CUMULATIVE_GCC[1] <- df$GCC[1]
  for (row in 2:nrow(df)){
    df$CUMULATIVE_GCC[row] <- df$GCC[row]+df$CUMULATIVE_GCC[row-1]
  }
  df <- cbind(df, CDF = 0)
  df$CDF <- df$CUMULATIVE_GCC / max(df$CUMULATIVE_GCC,na.rm=TRUE)
  df
}

# Initialize data frame
phenocam_cdf <- data.frame(matrix(ncol = 8, nrow = 0))
x <- c("phenocam_name", "eco_state", "year", "GCC", "DOY", "CUMULATIVE_GCC", "CDF", "baseline")
colnames(phenocam_cdf) <- x

# Calculate cumulative GCC and CDF for all phenocams and years
for (camera in unique(phenocam_data$phenocam_name)){
  for (year in unique(phenocam_data[phenocam_data$phenocam_name==camera,]$year)){
    subset <- phenocam_data[phenocam_data$phenocam_name==camera & phenocam_data$year==year,]
    eco_state <- subset$eco_state[1]
    #start_doy <- season_start_and_end[season_start_and_end$phenocam_name==camera & season_start_and_end$year==year,]$SOS_25
    #end_doy <- season_start_and_end[season_start_and_end$phenocam_name==camera & season_start_and_end$year==year,]$EOS_25
    start_doy <- min(subset$doy)
    end_doy <- max(subset$doy)
    baseline <- 0
    if (start_doy>=1 & start_doy<=366 & end_doy>=1 & end_doy<=366){
      subset <- cumulative_GCC(subset$gcc_scaled_yearly,subset$doy,start_doy,end_doy,baseline)
      subset <- cbind(phenocam_name = camera, eco_state = eco_state, year = year, subset, baseline = baseline)
      phenocam_cdf <- rbind(phenocam_cdf, subset)
    }
    #subset <- cumulative_GCC(subset$gcc_smooth,subset$doy,start_doy,end_doy,baseline)
    #subset <- cbind(phenocam_name = camera, year = year, subset, baseline = baseline)
    #phenocam_cdf <- rbind(phenocam_cdf, subset)
  }
}

rm(x,camera,year,subset,eco_state,start_doy,end_doy,baseline)

#################################
# Record CDF at 60-day intervals
#################################

phenocam_cdf_comparison <- data.frame(matrix(ncol = 9, nrow = 0))
x <- c("phenocam_name", "eco_state", "year",
       "CDF_30", "CDF_90", "CDF_150", "CDF_210", "CDF_270", "CDF_330")
colnames(phenocam_cdf_comparison) <- x

for (camera in unique(phenocam_cdf$phenocam_name)){
  for (year in unique(phenocam_cdf[phenocam_cdf$phenocam_name==camera,]$year)){
    subset <- phenocam_cdf[phenocam_cdf$phenocam_name==camera & phenocam_cdf$year==year,]
    if(subset$DOY[1] == 1){
      eco_state <- subset$eco_state[1]
      #rainfall <- subset$rainfall[1]
      cdf_30 <- subset[subset$DOY==30,]$CDF
      cdf_90 <- subset[subset$DOY==90,]$CDF
      cdf_150 <- subset[subset$DOY==150,]$CDF
      cdf_210 <- subset[subset$DOY==210,]$CDF
      cdf_270 <- subset[subset$DOY==270,]$CDF
      cdf_330 <- subset[subset$DOY==330,]$CDF
      new_row <- data.frame(phenocam_name = camera,
                            eco_state = eco_state,
                            year = year,
                            CDF_30 = cdf_30,
                            CDF_90 = cdf_90,
                            CDF_150 = cdf_150,
                            CDF_210 = cdf_210,
                            CDF_270 = cdf_270,
                            CDF_330 = cdf_330)
      phenocam_cdf_comparison <- rbind(phenocam_cdf_comparison, new_row)
    }
  }
}

rm(x,camera,year,subset,eco_state,cdf_30,cdf_90,
   cdf_210,cdf_270,cdf_330,new_row)

###############
# SAVE RESULTS
###############

save(list=c("phenocam_data","season_start_and_end","phenocam_cdf"), file="outputs/phenocam_processed.RData")