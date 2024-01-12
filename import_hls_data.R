## Code to import HLS data into R data frame, clean data, and calculate VIs

# Set working directory
setwd("~/hls-phenocam-plots")

# Load packages
#none so far

# Define functions

# Calculates the GCC
calc_GCC <- function(blue,green,red){
  GCC <- green/(blue+green+red)
  GCC
}

# Calculates the NDVI
calc_NDVI <- function(red,NIR){
  NDVI <- (NIR-red)/(NIR+red)
}

# Calculates the EVI
calc_EVI <- function(blue,red,NIR, L=1, C1=6, C2=7.5, G=2.5){
  EVI <- G*(NIR-red)/(NIR + C1*red - C2*blue + L)
}

# Calculates the value (0 or 1) of a bit at given position, given an HLS Quality
# number in decimal format
bit_value <- function(QA,bit_number){
  quotient <- QA%/%(2^bit_number)
  value <- quotient-((quotient%/%2)*2)
  value
}

# Calculates the date from year and day-of-year, returns in date format
calc_date <- function(year,doy){
  year <- as.character(year)
  date <- as.Date(doy-1,origin=paste(year,"-01-01",sep=""))
  date
}

# Returns ecological site and state based on phenocam name
return_eco_state <- function(phenocam_name){
  if(phenocam_name == "jershrubland"){
    eco_state = "Sandy shrubland"
  } else if (phenocam_name == "jershrubland2"){
    eco_state = "Sandy shrubland"
  } else if (phenocam_name == "jerbajada"){
    eco_state = "Gravelly shrubland"
  } else if (phenocam_name == "jernovel"){
    eco_state = "Sandy shrub-invaded grassland"
  } else if (phenocam_name == "jernovel2"){
    eco_state = "Sandy shrub-invaded grassland"
  } else if (phenocam_name == "jergrassland"){
    eco_state = "Sandy shrub-invaded grassland"
  } else if (phenocam_name == "jergrassland2"){
    eco_state = "Sandy shrub-invaded grassland"
  } else if (phenocam_name == "jernort"){
    eco_state = "Sandy shrubland"
  } else if (phenocam_name == "ibp"){
    eco_state = "Sandy shrub-invaded grassland"
  } else if (phenocam_name == "jernwern"){
    eco_state = "Sandy shrubland"
  } else if (phenocam_name == "NEON.D14.JORN.DP1.00033") {
    eco_state = "Sandy shrub-invaded grassland"
  } else if (phenocam_name == "jersand"){
    eco_state = "Gravelly shrubland"
  } else{
    eco_state = "unknown"
  }
  eco_state
}

# Imports a CSV of HLS data (produced by Python code import_HLS_pixel_data.py),
# filters out NA values, calculates vegetation indices, calculates QA bit
# values, and returns the outputs in a new data frame
hls_csv_to_df <- function(filename){
  hls_data <- read.csv(filename)
  hls_data <- hls_data[!is.na(hls_data$Red_mean),] # Remove NA values
  processed_data <- data.frame(DATE = calc_date(hls_data$Year, hls_data$DOY),
                               YEAR = hls_data$Year,
                               DOY = hls_data$DOY,
                               PHENOCAM_NAME = hls_data$Phenocam,
                               ECO_STATE = sapply(hls_data$Phenocam,return_eco_state),
                               SATELLITE = hls_data$Satellite,
                               GCC = calc_GCC(hls_data$Blue_mean, hls_data$Green_mean, hls_data$Red_mean),
                               NDVI = calc_NDVI(hls_data$Red_mean, hls_data$NIRNarrow_mean),
                               EVI = calc_EVI(hls_data$Blue_mean, hls_data$Red_mean, hls_data$NIRNarrow_mean),
                               STD_B = hls_data$Blue_std,
                               STD_G = hls_data$Green_std,
                               STD_R = hls_data$Red_std,
                               STD_NIR = hls_data$NIRNarrow_std,
                               Aerosol1 = bit_value(hls_data$Quality, 7),
                               Aerosol2 = bit_value(hls_data$Quality, 6),
                               Water = bit_value(hls_data$Quality, 5),
                               Snow = bit_value(hls_data$Quality, 4),
                               CloudShadow = bit_value(hls_data$Quality, 3),
                               CloudAdjacent = bit_value(hls_data$Quality, 2),
                               Cloud = bit_value(hls_data$Quality, 1))
  processed_data
}

# Put HLS data from all phenocams into one dataframe
hls_data_center <- rbind(hls_csv_to_df('data/outputs_HLS/jershrubland_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/jershrubland2_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/jernovel_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/jernovel2_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/jergrassland_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/jergrassland2_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/jerbajada_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/jernort_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/ibp_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/jernwern_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/NEON.D14.JORN.DP1.00033_center.csv'),
                         hls_csv_to_df('data/outputs_HLS/jersand_center.csv'))

hls_data_north <- rbind(hls_csv_to_df('data/outputs_HLS/jershrubland_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/jershrubland2_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/jernovel_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/jernovel2_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/jergrassland_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/jergrassland2_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/jerbajada_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/jernort_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/ibp_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/jernwern_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/NEON.D14.JORN.DP1.00033_north.csv'),
                         hls_csv_to_df('data/outputs_HLS/jersand_north.csv'))

# Create a subset of HLS data free from high aerosol, cloud, cloud shadow, and cloud adjacent
hls_center_clean <- hls_data_center[!(hls_data_center$Aerosol1==1 & hls_data_center$Aerosol2==1)
                                    & hls_data_center$Cloud==0
                                    & hls_data_center$CloudShadow==0
                                    & hls_data_center$CloudAdjacent==0,]

hls_north_clean <- hls_data_north[!(hls_data_north$Aerosol1==1 & hls_data_north$Aerosol2==1)
                                  & hls_data_north$Cloud==0
                                  & hls_data_north$CloudShadow==0
                                  & hls_data_north$CloudAdjacent==0,]

# Sort data (by date, and then sensor)
hls_center_clean <- hls_center_clean[order(hls_center_clean$DATE),]
hls_north_clean <- hls_north_clean[order(hls_north_clean$DATE),]

# Save results
save(list=c("hls_center_clean","hls_north_clean"), file="outputs/hls_v20_imported.RData")