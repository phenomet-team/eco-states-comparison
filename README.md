# Code from: Novel use of image time-series to distinguish dryland vegetation responses to wet and dry years

## Directories
- data: contains data products derived from PhenoCam, HLS, and precipitation raw data
- outputs: contains outputs of the following code in RData format

## General notes:
- In all scripts, please check directory names to ensure that the directories on your machine match the directories in the code

## Data import scripts
**download_phenocam_data.R** 
 - Downloads phenocam data and imports it into dataframe (output: phenocam_data_3day)

**import_HLS_pixel_data.py**
- Imports pixel values surrounding phenocams in study, outputs CSV and PKL files 

**import_hls_data.R**
- Imports HLS v2.0 data from csv files stored on the computer
- Calculates vegetation indices and quality bit values, outputs cleaned HLS dataframes (hls_center_clean, hls_north_clean)
- Saves output dataframes into RData file ("outputs/hls_v20_processed.RData")

## Data processing scripts (to be run after import)
**process_hls_data.R**
- Run only after running import_hls_v20_data.R at least once
- Performs LOESS smoothing and scaling on HLS v2.0 EVI, calculates season start and end dates, attaches rainfall data, calculates cumulative EVI
- Saves output dataframes into RData file (hls_smooth_scaled, season_start_and_end_hls, hls_cdf -> "outputs/hls_v14_processed.RData")

**process_phenocam_data.R**
- Run only after running download_phenocam_data.R
- Performs scaling on PhenoCam GCC, calculates season start and end dates, calculates cumulative GCC
- Saves output dataframes into RData file (phenocam_data, season_start_and_end, phenocam_cdf -> "outputs/phenocam_processed.RData")

**download_analyze_precip_data.R**
- Requires precipitation CSV files
- Aggregates precipitation (both standard and tipping bucket gauges) into rain-year totals

## Data analysis scripts (to be run after processing)
**eco_states_paper_analysis.R**
- Performs main statistical analysis tasks
- Generates figures for paper
