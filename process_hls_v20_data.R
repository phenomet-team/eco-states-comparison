## Code to process HLS time-series data

## NOTE: Before running this code, please run import_hls_data.R 
## or load the previously-imported HLS data

# Set working directory
setwd("~/hls-phenocam-plots")

# (Optional) Load previously-imported HLS data
load("outputs/hls_v20_imported.RData")

# Load packages
library(lubridate)

#############################
# SMOOTH THE HLS DATA (LOESS)
#############################

hls_data <- hls_center_clean[hls_center_clean$YEAR>=2016,]
hls_data <- cbind(hls_data,EVI_smooth=0)
hls_smooth <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(hls_smooth) <- c('DATE', 'YEAR', 'DOY', 'PHENOCAM_NAME', 'ECO_STATE',
                          'LOESS_SPAN', 'EVI_smooth')
loess_span <- 0.03

for (camera in unique(hls_data$PHENOCAM_NAME)){
  subset <- hls_data[hls_data$PHENOCAM_NAME == camera,]
  eco_state <- subset$ECO_STATE[1]
  subset$DATE <- as.numeric(subset$DATE)
  min_date <- min(subset$DATE, na.rm = TRUE)
  max_date <- max(subset$DATE, na.rm = TRUE)
  date_range <- seq(min_date,max_date,1)
  year <- year(as.Date(date_range, origin = "1970-01-01"))
  doy <- yday(as.Date(date_range, origin = "1970-01-01"))
  EVI_smooth <- predict(loess(EVI ~ DATE, data=subset, span=loess_span, control=loess.control(surface="direct")),newdata=date_range)
  subset_smooth <- data.frame(DATE = as.Date(date_range,origin="1970-01-01"),
                              YEAR = year,
                              DOY = doy,
                              PHENOCAM_NAME = camera,
                              ECO_STATE = eco_state,
                              LOESS_SPAN = loess_span,
                              EVI_smooth = EVI_smooth)
  hls_smooth <- rbind(hls_smooth,subset_smooth)
}

rm(camera,date_range,doy,EVI_smooth,loess_span,eco_state,max_date,min_date,year,subset,subset_smooth)

###############################################
# SCALE THE DATA (NORMALIZE BY SENSOR AND YEAR)
###############################################

hls_smooth_scaled <- cbind(hls_smooth,EVI_scaled_yearly = 0)

for (camera in unique(hls_smooth$PHENOCAM_NAME)){
  for (year in unique(hls_smooth[hls_smooth$PHENOCAM_NAME==camera,]$YEAR)){
    subset <- hls_smooth[hls_smooth$PHENOCAM_NAME==camera & hls_smooth$YEAR==year,]
    EVI_max <- max(subset$EVI_smooth, na.rm = TRUE)
    EVI_min <- min(subset$EVI_smooth, na.rm = TRUE)
    subset$EVI_scaled_yearly <- (subset$EVI_smooth - EVI_min) / (EVI_max - EVI_min)
    hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==camera & hls_smooth_scaled$YEAR==year,]$EVI_scaled_yearly <- subset$EVI_scaled_yearly
  }
}

rm(camera,EVI_max,EVI_min,year,subset)

#######################################
# Calculate season start and end dates
#######################################

# NOTES: For thresholding, added the restrictions that SOS must occur after GCC
# reaches its pre-peak minimum value and EOS must occur before GCC reaches its
# post-peak minimum value

# Initialize data frame
season_start_and_end_hls <- data.frame(matrix(ncol = 15, nrow = 0))
x <- c("phenocam_name", "eco_state", "year", "peak_EVI",
       "DOY_peak", "DOY_min_pre_peak", "DOY_min_post_peak",
       "SOS_10", "SOS_15", "SOS_25", "SOS_50", "EOS_10", "EOS_15", "EOS_25", "EOS_50")
colnames(season_start_and_end_hls) <- x

for (camera in unique(hls_smooth_scaled$PHENOCAM_NAME)){
  for (year in unique(hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==camera,]$YEAR)){
    subset <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==camera & hls_smooth_scaled$YEAR==year,]
    eco_state <- subset$ECO_STATE[1]
    # Peak - day of the year when GCC is at its maximum value
    peak <- min(subset[subset$EVI_scaled_yearly == 1,]$DOY, na.rm = TRUE)
    peak_EVI <- subset[subset$DOY==peak,]$EVI_smooth
    # Pre-peak minimum - day of the year when GCC is at its minimum value, prior to peak
    min1 <- min(subset[subset$DOY < peak,]$EVI_scaled_yearly, na.rm = TRUE)
    min1 <- min(subset[subset$EVI_scaled_yearly == min1,]$DOY, na.rm = TRUE)
    # Post-peak minimum - day of the year when GCC is at its minimum value, post peak
    min2 <- min(subset[subset$DOY > peak,]$EVI_scaled_yearly, na.rm = TRUE)
    min2 <- max(subset[subset$EVI_scaled_yearly == min2,]$DOY, na.rm = TRUE)
    # Start of season - first day of the year when GCC is above the threshold
    SOS_10 <- min(subset[subset$EVI_scaled_yearly > 0.10 & subset$DOY >= min1,]$DOY, na.rm = TRUE)
    SOS_15 <- min(subset[subset$EVI_scaled_yearly > 0.15 & subset$DOY >= min1,]$DOY, na.rm = TRUE)
    SOS_25 <- min(subset[subset$EVI_scaled_yearly > 0.25 & subset$DOY >= min1,]$DOY, na.rm = TRUE)
    SOS_50 <- min(subset[subset$EVI_scaled_yearly > 0.50 & subset$DOY >= min1,]$DOY, na.rm = TRUE)
    # End of season - last day of the year when GCC is above the threshold
    EOS_10 <- max(subset[subset$EVI_scaled_yearly > 0.10 & subset$DOY <= min2,]$DOY, na.rm = TRUE)
    EOS_15 <- max(subset[subset$EVI_scaled_yearly > 0.15 & subset$DOY <= min2,]$DOY, na.rm = TRUE)
    EOS_25 <- max(subset[subset$EVI_scaled_yearly > 0.25 & subset$DOY <= min2,]$DOY, na.rm = TRUE)
    EOS_50 <- max(subset[subset$EVI_scaled_yearly > 0.50 & subset$DOY <= min2,]$DOY, na.rm = TRUE)
    new_row <- data.frame(phenocam_name = camera,
                          eco_state = eco_state,
                          year = year,
                          peak_EVI = peak_EVI,
                          DOY_peak = peak,
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
    season_start_and_end_hls <- rbind(season_start_and_end_hls,new_row)
  }
}

# Remove 2013-2015 (not enough S30)
season_start_and_end_hls <- season_start_and_end_hls[season_start_and_end_hls$year>=2016 & season_start_and_end_hls$year<2023,]

rm(x,camera,year,subset,eco_state,peak_EVI,peak,min1,min2,
   SOS_10,SOS_15,SOS_25,SOS_50,EOS_10,EOS_15,EOS_25,EOS_50,new_row)

###########################
# Calculate cumulative EVI
###########################

# Define function
cumulative_EVI <- function(evi_df,time_df,start_doy,end_doy,baseline=0){
  df <- data.frame(EVI = evi_df, DOY = time_df)
  df <- cbind(df, CUMULATIVE_EVI = 0)
  df <- df[df$DOY >= start_doy & df$DOY <= end_doy,]
  df$EVI <- df$EVI - baseline
  # commented-out line below does not work as intended, but is supposed to make
  # it so that if GCC dips below baseline, cumulative GCC does not decrease
  # df[df$GCC < 0,]$GCC <- 0
  df$CUMULATIVE_EVI[1] <- df$EVI[1]
  for (row in 2:nrow(df)){
    df$CUMULATIVE_EVI[row] <- df$EVI[row]+df$CUMULATIVE_EVI[row-1]
  }
  df <- cbind(df, CDF = 0)
  df$CDF <- df$CUMULATIVE_EVI / max(df$CUMULATIVE_EVI,na.rm=TRUE)
  df
}

# Initialize data frame
hls_cdf <- data.frame(matrix(ncol = 8, nrow = 0))
x <- c("phenocam_name", "eco_state", "year", "EVI", "DOY", "CUMULATIVE_EVI", "CDF", "baseline")
colnames(hls_cdf) <- x

# Calculate cumulative EVI and CDF for all locations and years
for (camera in unique(hls_smooth_scaled$PHENOCAM_NAME)){
  for (year in unique(season_start_and_end_hls[season_start_and_end_hls$phenocam_name==camera,]$year)){
    subset <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==camera & hls_smooth_scaled$YEAR==year,]
    eco_state <- subset$ECO_STATE[1]
    #start_doy <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==camera & season_start_and_end_hls$year==year,]$SOS_25
    #end_doy <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==camera & season_start_and_end_hls$year==year,]$EOS_25
    start_doy <- min(subset$DOY)
    end_doy <- max(subset$DOY)
    baseline <- 0
    if (start_doy>=1 & start_doy<=366 & end_doy>=1 & end_doy<=366){
      subset <- cumulative_EVI(subset$EVI_smooth,subset$DOY,start_doy,end_doy,baseline)
      subset <- cbind(data.frame(phenocam_name = camera,
                                 eco_state = eco_state,
                                 year = year),
                      subset,
                      data.frame(baseline = baseline))
      hls_cdf <- rbind(hls_cdf, subset)
    }
  }
}

rm(x,camera,year,subset,eco_state,start_doy,end_doy,baseline)

############################################
# Record cumulative EVI at 60-day intervals
############################################

hls_cdf_comparison <- data.frame(matrix(ncol = 10, nrow = 0))
x <- c("phenocam_name", "eco_state", "year", "rainfall",
       "CDF_30", "CDF_90", "CDF_150", "CDF_210", "CDF_270", "CDF_330")
colnames(hls_cdf_comparison) <- x

for (camera in unique(hls_cdf$phenocam_name)){
  for (year in unique(hls_cdf[hls_cdf$phenocam_name==camera,]$year)){
    subset <- hls_cdf[hls_cdf$phenocam_name==camera & hls_cdf$year==year,]
    eco_state <- subset$eco_state[1]
    cum_EVI_30 <- subset[subset$DOY==30,]$CUMULATIVE_EVI
    cum_EVI_90 <- subset[subset$DOY==90,]$CUMULATIVE_EVI
    cum_EVI_150 <- subset[subset$DOY==150,]$CUMULATIVE_EVI
    cum_EVI_210 <- subset[subset$DOY==210,]$CUMULATIVE_EVI
    cum_EVI_270 <- subset[subset$DOY==270,]$CUMULATIVE_EVI
    cum_EVI_330 <- subset[subset$DOY==330,]$CUMULATIVE_EVI
    new_row <- data.frame(phenocam_name = camera,
                          eco_state = eco_state,
                          year = year,
                          CDF_30 = cum_EVI_30,
                          CDF_90 = cum_EVI_90,
                          CDF_150 = cum_EVI_150,
                          CDF_210 = cum_EVI_210,
                          CDF_270 = cum_EVI_270,
                          CDF_330 = cum_EVI_330)
    hls_cdf_comparison <- rbind(hls_cdf_comparison, new_row)
  }
}

rm(x,camera,year,subset,eco_state,cum_EVI_30,cum_EVI_90,cum_EVI_150,
   cum_EVI_210,cum_EVI_270,cum_EVI_330,new_row)

###############
# SAVE RESULTS
###############

save(list=c("hls_smooth_scaled","season_start_and_end_hls","hls_cdf","hls_cdf_comparison"), file="outputs/hls_v20_processed.RData")