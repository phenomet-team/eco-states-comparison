## Code to download Phenocam time-series data and import it into data frames

# Set directory
setwd("~/hls-phenocam-plots")

# Load packages
library(phenocamr)

# Define functions
phenocam_csv_to_df <- function(list, interval, Skip=24){
  if(interval == 1){
    filename <- list[2]
  }
  else{
    filename <- list[1]
  }
  PhenocamName <- list[3]
  EcoState <- list[4]
  phenocam_data <- read.csv(filename, skip=Skip)
  phenocam_data$date <- as.Date(phenocam_data$date, "%Y-%m-%d")
  min_date <- min(phenocam_data[complete.cases(phenocam_data$gcc_90),]$date)
  max_date <- max(phenocam_data[complete.cases(phenocam_data$gcc_90),]$date)
  phenocam_data <- phenocam_data[phenocam_data$date>=min_date & phenocam_data$date<=max_date,]
  #phenocam_data <- phenocam_data[complete.cases(phenocam_data$gcc_90),]
  phenocam_data <- cbind(phenocam_data,eco_state=EcoState)
  phenocam_data <- cbind(phenocam_data,phenocam_name=PhenocamName)
  phenocam_data
}

# Download 3-day phenocam GCC
freq <- "3"
download_phenocam(site = 'jershrubland$', veg_type = 'XX', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'jershrubland2$', veg_type='XX', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'jernovel$', veg_type='XX', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'jernovel2$', veg_type = 'XX', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'jergrassland$', veg_type = 'XX', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'jergrassland2$', veg_type='XX', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'jerbajada$', veg_type = 'SH', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'jernort$', veg_type = 'XX', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'ibp$', veg_type = 'XX', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'jernwern$', veg_type = 'XX', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'NEON.D14.JORN.DP1.00033$', veg_type = 'GR', frequency = freq, out_dir = './data/phenocam/')
download_phenocam(site = 'jersand$', veg_type = 'SH', frequency = freq, out_dir = './data/phenocam/')

# Site information
# 3-day filename, 1-day filename, site name, site code
jershrubland  <- c('./data/phenocam/jershrubland_XX_1000_3day.csv',
                   './data/phenocam/jershrubland_XX_1000_1day.csv',
                   'jershrubland', "Sandy shrubland")
jershrubland2 <- c('./data/phenocam/jershrubland2_XX_1000_3day.csv',
                   './data/phenocam/jershrubland2_XX_1000_1day.csv',
                   'jershrubland2', "Sandy shrubland")
jernovel      <- c('./data/phenocam/jernovel_XX_1000_3day.csv',
                   './data/phenocam/jernovel_XX_1000_1day.csv',
                   'jernovel', "Sandy shrub-invaded grassland")
jernovel2     <- c('./data/phenocam/jernovel2_XX_1000_3day.csv',
                   './data/phenocam/jernovel2_XX_1000_1day.csv',
                   'jernovel2', "Sandy shrub-invaded grassland")
jergrassland  <- c('./data/phenocam/jergrassland_XX_1000_3day.csv',
                   './data/phenocam/jergrassland_XX_1000_1day.csv',
                   'jergrassland', "Sandy shrub-invaded grassland")
jergrassland2 <- c('./data/phenocam/jergrassland2_XX_1000_3day.csv',
                   './data/phenocam/jergrassland2_XX_1000_1day.csv',
                   'jergrassland2', "Sandy shrub-invaded grassland")
jerbajada     <- c('./data/phenocam/jerbajada_SH_1000_3day.csv',
                   './data/phenocam/jerbajada_SH_1000_1day.csv',
                   'jerbajada', "Gravelly shrubland")
jernort       <- c('./data/phenocam/jernort_XX_1000_3day.csv',
                   './data/phenocam/jernort_XX_1000_1day.csv',
                   'jernort', "Sandy shrubland")
ibp           <- c('./data/phenocam/ibp_XX_1000_3day.csv',
                   './data/phenocam/ibp_XX_1000_1day.csv',
                   'ibp', "Sandy shrub-invaded grassland")
jernwern      <- c('./data/phenocam/jernwern_XX_1000_3day.csv',
                   './data/phenocam/jernwern_XX_1000_1day.csv',
                   'jernwern', "Sandy shrubland")
neon          <- c('./data/phenocam/NEON.D14.JORN.DP1.00033_GR_1000_3day.csv',
                   './data/phenocam/NEON.D14.JORN.DP1.00033_GR_1000_1day.csv',
                   'NEON.D14.JORN.DP1.00033', "Sandy shrub-invaded grassland")
jersand       <- c('./data/phenocam/jersand_SH_1000_3day.csv',
                   './data/phenocam/jersand_SH_1000_1day.csv',
                   'jersand', "Gravelly shrubland")

phenocam_data_3day <- rbind(phenocam_csv_to_df(jershrubland,interval=3),
                            phenocam_csv_to_df(jershrubland2,interval=3),
                            phenocam_csv_to_df(jernovel,interval=3),
                            phenocam_csv_to_df(jernovel2,interval=3),
                            phenocam_csv_to_df(jergrassland,interval=3),
                            phenocam_csv_to_df(jergrassland2,interval=3),
                            phenocam_csv_to_df(jerbajada,interval=3),
                            phenocam_csv_to_df(jernort,interval=3),
                            phenocam_csv_to_df(ibp,interval=3),
                            phenocam_csv_to_df(jernwern,interval=3),
                            phenocam_csv_to_df(neon,interval=3),
                            phenocam_csv_to_df(jersand,interval=3))

# Remote extraneous variables from the workspace
rm(freq, jershrubland, jershrubland2, jernovel, jernovel2, jergrassland,
   jergrassland2, jerbajada, jernort, ibp, jernwern, neon, jersand)

# Save output
save("phenocam_data_3day", file="outputs/phenocam_imported.RData")