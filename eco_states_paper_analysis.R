
# Set working directory
setwd("~/hls-phenocam-plots")

# Load packages
library(dplyr)
detach("package:dplyr", unload = TRUE)
library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)
library(emmeans)
library(multcomp)

# Load data (PhenoCam, HLS v1.4, and rainfall)
load("outputs/hls_v20_processed.RData")
load("outputs/phenocam_processed.RData")
load("outputs/annual_precip.RData")

# Combine precip data with season start and end data
season_start_and_end <- left_join(season_start_and_end,phenocam_rain_gauge_pairs,by = "phenocam_name")
season_start_and_end <- left_join(season_start_and_end,select(annual_site_precip_selected_years,c(site,year,wet_dry_ave_spring_summer,wet_dry_ave_winter)), by = c("site","year"))
season_start_and_end <- season_start_and_end %>% rename(rainfall = wet_dry_ave_spring_summer)
season_start_and_end_hls <- left_join(season_start_and_end_hls,phenocam_rain_gauge_pairs,by = "phenocam_name")
season_start_and_end_hls <- left_join(season_start_and_end_hls,select(annual_site_precip_selected_years,c(site,year,wet_dry_ave_spring_summer,wet_dry_ave_winter)), by = c("site","year"))
season_start_and_end_hls <- season_start_and_end_hls %>% rename(rainfall = wet_dry_ave_spring_summer)

# Remove unwanted data
# 2013
season_start_and_end <- season_start_and_end[season_start_and_end$year > 2013,]
# Inf values (jernwern 2017)
season_start_and_end <- season_start_and_end[season_start_and_end$SOS_25 != Inf & season_start_and_end$EOS_25 != Inf,]

#########################
# FIGURE: PRECIPITATION
#########################

# Total spring-summer precipitation for all relevant gauges
ggplot(data = annual_site_precip_selected_years, aes(x=year, y=total_rainfall_spring_summer, color=site)) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = 203.3, linetype = "Wet threshold")) +
  geom_hline(aes(yintercept = 139.1, linetype = "Dry threshold")) +
  scale_linetype_manual(values = c("dashed","dotdash")) +
  xlab("Calendar year") +
  ylab("Total spring and summer precipitation (mm)") +
  labs(color="PhenoCams", linetype=" ") +
  scale_color_discrete(labels=c(' \njernovel\njernovel2\njershrubland\njershrubland2\n ', 'jergrassland2', 
                               'jernwern', 'jernort', 'jerbajada', 'jersand', ' \njergrassland\nibp\nNEON\n ')) +
  theme_bw()
ggsave("figures/precipitation_plot.jpg", height = 8, width = 10, units = "in", dpi = 300)

###########################################
# FIGURE: Phenocam-HLS SOS-EOS comparisons
###########################################

# Merge data for sensor comparison scatterplots (25% threshold only)
sos_eos_phen <- season_start_and_end %>% select(phenocam_name,eco_state,year,rainfall,peak,SOS_25,EOS_25)
sos_eos_phen <- sos_eos_phen %>% rename(peak_phen = peak, SOS_phen = SOS_25, EOS_phen = EOS_25)
sos_eos_hls <- season_start_and_end_hls %>% select(phenocam_name,eco_state,year,rainfall,peak_EVI,DOY_peak,SOS_25,EOS_25)
sos_eos_hls <- sos_eos_hls %>% rename(peak_hls = DOY_peak, SOS_hls = SOS_25, EOS_hls = EOS_25)
sos_eos_merge <- inner_join(sos_eos_phen, sos_eos_hls, by=c('phenocam_name','year','eco_state','rainfall'))
sos_eos_merge <- sos_eos_merge %>% mutate(rainfall_factor = factor(rainfall, levels = c("dry", "average", "wet")))
rm(sos_eos_phen,sos_eos_hls)

# Generate scatter plots

p1 <- ggplot(data = sos_eos_merge, aes(x = SOS_phen, y = SOS_hls, shape = eco_state, color = rainfall_factor)) +
  geom_point(size = 3) +
  geom_abline(intercept=0,slope=1) +
  xlim(0,250) +
  ylim(0,250) +
  xlab("Phenocam SOS - 25% threshold") +
  ylab("HLS SOS - 25% threshold") +
  labs(shape="Ecological State", color="Rainfall") +
  theme_bw() +
  scale_color_discrete(breaks=c('wet','average','dry')) +
  ggtitle("(a) SOS")

p3 <- ggplot(data = sos_eos_merge, aes(x = EOS_phen, y = EOS_hls, shape = eco_state, color = rainfall_factor)) +
  geom_point(size = 3) +
  geom_abline(intercept=0,slope=1) +
  xlim(150,400) +
  ylim(150,400) +
  xlab("Phenocam EOS - 25% threshold") +
  ylab("HLS EOS - 25% threshold") +
  ggtitle("(c) EOS") +
  theme_bw() +
  scale_color_discrete(breaks=c('wet','average','dry')) +
  theme(legend.position='none')

p2 <- ggplot(data = sos_eos_merge, aes(x = peak_phen, y = peak_hls, shape = eco_state, color = rainfall_factor)) +
  geom_point(size = 3) +
  geom_abline(intercept=0,slope=1) +
  xlim(75,325) +
  ylim(75,325) +
  xlab("Phenocam Peak DOY - 25% threshold") +
  ylab("HLS Peak DOY - 25% threshold") +
  ggtitle("(b) Peak DOY") +
  theme_bw() +
  scale_color_discrete(breaks=c('wet','average','dry')) +
  theme(legend.position='none')

p4 <- get_legend(p1)

p1 <- p1 + theme(legend.position='none')

plot_grid(p1,p2,p3,p4, nrow=2, ncol=2)

ggsave("figures/sos_peak_eos_scatter.jpg", height = 8.25, width = 7.91, units = "in", dpi = 300, bg="white")

######################################################
# SUPPLEMENTAL FIGURES: VISUAL INSPECTION OF OUTLIERS
######################################################

# SOS

site <- "jerbajada"
number <- "11"
year <- 2022
label <- "(1)"
transition <- "SOS"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_SOS <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$SOS_25
phen_SOS <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$SOS_25
eco_state <- subset1$ECO_STATE[1]
a_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_SOS, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
a_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_SOS, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "jersand"
number <- "12"
year <- 2016
label <- "(2)"
transition <- "SOS"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_SOS <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$SOS_25
phen_SOS <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$SOS_25
eco_state <- subset1$ECO_STATE[1]
b_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_SOS, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
b_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_SOS, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "jersand"
number <- "12"
year <- 2017
label <- "(3)"
transition <- "SOS"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_SOS <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$SOS_25
phen_SOS <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$SOS_25
eco_state <- subset1$ECO_STATE[1]
c_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_SOS, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
c_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_SOS, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "jersand"
number <- 12
year <- 2018
label <- "(4)"
transition <- "SOS"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_SOS <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$SOS_25
phen_SOS <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$SOS_25
eco_state <- subset1$ECO_STATE[1]
d_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_SOS, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
d_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_SOS, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "NEON.D14.JORN.DP1.00033"
number <- "6"
year <- 2017
label <- "(5)"
transition <- "SOS"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_SOS <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$SOS_25
phen_SOS <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$SOS_25
eco_state <- subset1$ECO_STATE[1]
e_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_SOS, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
e_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_SOS, color = "black", linetype = "dashed") +
  xlim(0,365) +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "NEON.D14.JORN.DP1.00033"
number <- 6
year <- 2018
label <- "(6)"
transition <- "SOS"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_SOS <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$SOS_25
phen_SOS <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$SOS_25
eco_state <- subset1$ECO_STATE[1]
f_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_SOS, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
f_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_SOS, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "NEON.D14.JORN.DP1.00033"
number <- "6"
year <- 2019
label <- "(7)"
transition <- "SOS"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_SOS <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$SOS_25
phen_SOS <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$SOS_25
eco_state <- subset1$ECO_STATE[1]
g_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_SOS, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
g_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_SOS, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

h_1 <- ggplot() + theme_nothing()
h_2 <- ggplot() + theme_nothing()

a_1 + b_1 + a_2 + b_2 + c_1 + d_1 + c_2 + d_2 +
e_1 + f_1 + e_2 + f_2 + g_1 + h_1 + g_2 + h_2 + plot_layout(ncol=2)

ggsave("figures/outliers_SOS.jpg", height = 10, width = 9, units = "in", dpi = 300, bg="white")

# Peak

site <- "jerbajada"
number <- "11"
year <- 2016
label <- "(8)"
transition <- "Peak"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_peak <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$DOY_peak
phen_peak <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$peak
eco_state <- subset1$ECO_STATE[1]
a_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_peak, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
a_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_peak, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "jernwern"
number <- "8"
year <- 2018
label <- "(9)"
transition <- "Peak"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_peak <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$DOY_peak
phen_peak <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$peak
eco_state <- subset1$ECO_STATE[1]
b_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_peak, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
b_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_peak, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "jergrassland"
number <- "2"
year <- 2019
label <- "(10)"
transition <- "Peak"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_peak <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$DOY_peak
phen_peak <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$peak
eco_state <- subset1$ECO_STATE[1]
c_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_peak, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
c_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_peak, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "jernort"
number <- "7"
year <- 2019
label <- "(11)"
transition <- "Peak"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_peak <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$DOY_peak
phen_peak <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$peak
eco_state <- subset1$ECO_STATE[1]
d_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_peak, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
d_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_peak, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "jersand"
number <- "12"
year <- 2020
label <- "(12)"
transition <- "Peak"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_peak <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$DOY_peak
phen_peak <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$peak
eco_state <- subset1$ECO_STATE[1]
e_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_peak, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
e_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_peak, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "jernovel"
number <- "4"
year <- 2020
label <- "(13)"
transition <- "Peak"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_peak <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$DOY_peak
phen_peak <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$peak
eco_state <- subset1$ECO_STATE[1]
f_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_peak, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
f_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_peak, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

site <- "jernort"
number <- "7"
year <- 2022
label <- "(14)"
transition <- "Peak"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_peak <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$DOY_peak
phen_peak <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$peak
eco_state <- subset1$ECO_STATE[1]
g_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_peak, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
g_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_peak, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

a_1 + b_1 + a_2 + b_2 + c_1 + d_1 + c_2 + d_2 +
e_1 + f_1 + e_2 + f_2 + g_1 + h_1 + g_2 + h_2 + plot_layout(ncol=2)

ggsave("figures/outliers_peak.jpg", height = 10, width = 9, units = "in", dpi = 300, bg="white")

# EOS (only one site)

site <- "NEON.D14.JORN.DP1.00033"
number <- "6"
year <- 2020
label <- "(15)"
transition <- "EOS"

subset1 <- hls_smooth_scaled[hls_smooth_scaled$PHENOCAM_NAME==site & hls_smooth_scaled$YEAR==year,]
subset2 <- phenocam_data[phenocam_data$phenocam_name==site & phenocam_data$year==year,]
hls_EOS <- season_start_and_end_hls[season_start_and_end_hls$phenocam_name==site & season_start_and_end_hls$year==year,]$EOS_25
phen_EOS <- season_start_and_end[season_start_and_end$phenocam_name==site & season_start_and_end$year==year,]$EOS_25
eco_state <- subset1$ECO_STATE[1]
a_1 <- ggplot(data = subset1, aes(x=DOY,y=EVI_smooth)) +
  geom_line() +
  ylab("EVI") +
  geom_vline(xintercept=hls_EOS, color = "black", linetype = "dashed") +
  labs(x = NULL) +
  theme_bw() +
  theme(plot.title = element_text(size=11)) +
  ggtitle(paste(label," ",eco_state,", ",year," ",transition, sep=""))
a_2 <- ggplot(data = subset2, aes(x=doy,y=smooth_gcc_90)) +
  geom_line() +
  geom_vline(xintercept=phen_EOS, color = "black", linetype = "dashed") +
  ylab("GCC") +
  xlab("DOY") +
  theme_bw()

a_1 / a_2

ggsave("figures/outliers_EOS.jpg", height = 3, width = 5, units = "in", dpi = 300, bg="white")

#################################################
# PREPARE DATA FOR ANOVA comparison (25% threshold)
#################################################
sos_eos_phen <- season_start_and_end %>% dplyr::select(phenocam_name,eco_state,year,rainfall,peak,SOS_25,EOS_25)
sos_eos_phen <- cbind(sos_eos_phen, sensor="PhenoCam")
sos_eos_hls <- season_start_and_end_hls %>% dplyr::select(phenocam_name,eco_state,year,rainfall,DOY_peak,SOS_25,EOS_25)
sos_eos_hls <- sos_eos_hls %>% rename(peak = DOY_peak)
sos_eos_hls <- cbind(sos_eos_hls, sensor="HLS")
sos_eos_combine <- rbind(sos_eos_phen,sos_eos_hls)
rm(sos_eos_phen,sos_eos_hls)

#############
# RUN ANOVAS
#############

# Suppress scientific notation
options(scipen = 9999)
options(contrasts = c("contr.sum", "contr.poly")) # Necessary for conducting Type III sums of squares tests
options("contrasts")

write_csv(sos_eos_combine, "sos_eos_combine.csv")
write_csv(season_start_and_end_hls, "season_start_and_end_hls.csv")

# Calculate season duration
season.duration <- sos_eos_combine %>%
  mutate(duration = EOS_25 - SOS_25)

write_csv(season.duration, "season_duration.csv")

# Run ANOVA models
aov_sos_3state <- lm(SOS_25 ~ sensor*eco_state*rainfall, data = sos_eos_combine)
aov_peak_3state <- lm(peak ~ sensor*eco_state*rainfall, data = sos_eos_combine)
aov_eos_3state <- lm(EOS_25 ~ sensor*eco_state*rainfall, data = sos_eos_combine)
aov_duration_3state <- lm(duration ~ sensor*eco_state*rainfall, data = season.duration) # This is new
aov_peak_EVI_3state <- lm(peak_EVI ~ eco_state*rainfall, data = season_start_and_end_hls)

# Type III tests of fixed effects for each model
sos.type3 <- car::Anova(aov_sos_3state, type = "III") %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "model_effect") %>%
  mutate(response = "Start of Season") 

peak.type3 <- car::Anova(aov_peak_3state, type = "III") %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "model_effect") %>%
  mutate(response = "Peak") 

eos.type3 <- car::Anova(aov_eos_3state, type = "III") %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "model_effect") %>%
  mutate(response = "End of Season") 

duration.type3 <- car::Anova(aov_duration_3state, type = "III") %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "model_effect") %>%
  mutate(response = "Season Duration") 

evi.type3 <- car::Anova(aov_peak_EVI_3state, type = "III") %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "model_effect") %>%
  mutate(response = "Peak EVI") 

# Stack the Type III test results together
Type3Tests <- bind_rows(sos.type3, peak.type3) %>%
  bind_rows(eos.type3) %>%
  bind_rows(duration.type3) %>%
  bind_rows(evi.type3) %>%
  dplyr::filter(!(model_effect %in% c("(Intercept)", "Residuals"))) %>%
  dplyr::select(response, model_effect, Sum_Sq = `Sum Sq`, Df, F_stat = `F value`, P = `Pr(>F)`) %>%
  mutate(p_value = if_else(P < 0.0001, "p < 0.0001", paste0("p = ", format(round(P, digits = 4), nsmall = 4))),
         singificant = if_else(P > 0.05 & P <= 0.1, "marginal",
                               if_else(P <= 0.05, "significant", "")))

write_csv(Type3Tests, "Phenology ANOVA Type3 Tests of Fixed Effects.csv")

unique(Type3Tests$model_effect)

# Loop through the 3-way models and calculate least squares means

anova.results <- list(aov_sos_3state, aov_peak_3state, aov_eos_3state, aov_duration_3state)
anova.names <- c("Start of Season", "Peak", "End of Season", "Season Duration")

sos.peak.eos.LSMeans <- data.frame() # Empty data.frame to collect the results

for(i in 1:length(anova.results)) {
  
  # Main effects
  sensor.main.lsm <- emmeans(anova.results[[i]], ~ sensor, adjust = "Sidak") %>%
    multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
    mutate(model_effect = "sensor",
           effect = "sensor main effect")
  
  eco_state.main.lsm <- emmeans(anova.results[[i]], ~ eco_state, adjust = "Sidak") %>%
    multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
    mutate(model_effect = "eco_state",
           effect = "eco_state main effect")
  
  rainfall.main.lsm <- emmeans(anova.results[[i]], ~ rainfall, adjust = "Sidak") %>%
    multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
    mutate(model_effect = "rainfall",
           effect = "rainfall main effect")
  
  # Two-way interactions
  sensor.eco_state.lsm <- emmeans(anova.results[[i]], ~ sensor:eco_state, adjust = "Sidak") %>%
    multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
    mutate(model_effect = "sensor:eco_state",
           effect = "sensor * eco_state interaction")
  
  sensor.rainfall.lsm <- emmeans(anova.results[[i]], ~ sensor:rainfall, adjust = "Sidak") %>%
    multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
    mutate(model_effect = "sensor:rainfall",
           effect = "sensor * rainfall interaction")
  
  eco_state.rainfall.lsm <- emmeans(anova.results[[i]], ~ eco_state:rainfall, adjust = "Sidak") %>%
    multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
    mutate(model_effect = "eco_state:rainfall",
           effect = "eco_state * rainfall interaction")
  
  # Three-way interaction
  sensor.eco_state.rainfall.lsm <- emmeans(anova.results[[i]], ~ sensor:eco_state:rainfall, adjust = "Sidak") %>%
    multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
    mutate(model_effect = "sensor:eco_state:rainfall",
           effect = "sensor * eco_state * rainfall interaction")
  
  lsms.3way.models <- bind_rows(sensor.eco_state.rainfall.lsm, eco_state.rainfall.lsm) %>%
    bind_rows(sensor.rainfall.lsm) %>%
    bind_rows(sensor.eco_state.lsm) %>%
    bind_rows(rainfall.main.lsm) %>%
    bind_rows(eco_state.main.lsm) %>%
    bind_rows(sensor.main.lsm) %>%
    mutate(mean_SE = paste(format(round(emmean, digits = 1), nsmall = 1), "+/-", format(round(SE, digits = 1), nsmall = 1)) %>% trimws()) %>%
    mutate(response = anova.names[i]) %>%
    mutate(.group = trimws(.group)) %>%
    dplyr::select(response, model_effect, effect, sensor, eco_state, rainfall, emmean, SE, mean_SE, .group, df, lower.CL, upper.CL)
  
  sos.peak.eos.LSMeans <- rbind(sos.peak.eos.LSMeans, lsms.3way.models)
  
}

# Add in peak_EVI results; this is from HLS only (no sensor effect)
# Main effects
hls.eco_state.lsm <- emmeans(aov_peak_EVI_3state, ~ eco_state, adjust = "Sidak") %>%
  multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
  mutate(model_effect = "eco_state",
         effect = "eco_state main effect")

hls.rainfall.lsm <- emmeans(aov_peak_EVI_3state, ~ rainfall, adjust = "Sidak") %>%
  multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
  mutate(model_effect = "rainfall",
         effect = "rainfall main effect")

# Two-way interaction
hls.eco_state.rainfall.lsm <- emmeans(aov_peak_EVI_3state, ~ eco_state:rainfall, adjust = "Sidak") %>%
  multcomp::cld(Letters = LETTERS) %>% as_tibble() %>%
  mutate(model_effect = "eco_state:rainfall",
         effect = "eco_state * rainfall interaction")

# Stack Peak EVI LSMeans together
hls.lsms <- bind_rows(hls.eco_state.rainfall.lsm, hls.rainfall.lsm) %>%
  bind_rows(hls.eco_state.lsm) %>%
  mutate(mean_SE = paste(format(round(emmean, digits = 1), nsmall = 1), "+/-", format(round(SE, digits = 1), nsmall = 1)) %>% trimws()) %>%
  mutate(response = "Peak EVI") %>%
  mutate(.group = trimws(.group)) %>%
  dplyr::select(response, model_effect, effect, eco_state, rainfall, emmean, SE, mean_SE, .group, df, lower.CL, upper.CL)

# Combine all the LSMeans from all 5 ANOVA models
lsms.all <- bind_rows(sos.peak.eos.LSMeans, hls.lsms) %>%
  left_join(Type3Tests %>% dplyr::select(response, model_effect, p_value)) %>%
  mutate(effect_name = paste0(effect, " (", p_value, ")"))

# Export as CSV
write_csv(lsms.all %>% dplyr::select(-effect_name), "LSMeans from Phenology ANOVA models.csv", na = "")

# Calculate number of observations for each eco state/sensor/rainfall combo
number_of_observations <- ddply(sos_eos_combine, c("eco_state","rainfall","sensor"), summarize,
                                   eco_state = eco_state,
                                   rainfall = rainfall,
                                   sensor = sensor,
                                   number_of_samples = n())
number_of_observations <- unique(number_of_observations)

##########################################################
# FIGURE: SOS, PEAK, EOS, and DURATION LEAST SQUARE MEANS
##########################################################

sos_peak_eos_dur_lsmeans <- sos.peak.eos.LSMeans[sos.peak.eos.LSMeans$model_effect=="sensor:eco_state:rainfall",] %>%
  mutate(sensor_factor = factor(sensor, levels = c("PhenoCam", "HLS"))) %>%
  mutate(rainfall_factor = factor(rainfall, levels = c("dry", "average", "wet"))) %>%
  mutate(response_factor = factor(response, levels = c("Start of Season", "Peak", "End of Season", "Season Duration")))

plot = ggplot(sos_peak_eos_dur_lsmeans, aes(x = eco_state, y = emmean, col = rainfall_factor)) +
  theme_bw() +
  geom_point(position = position_dodge(0.9)) +
  facet_grid(sensor_factor ~ response_factor, scales = "free_x", labeller = labeller(response_factor = c(`Start of Season` = "(a) Start of Season", `Peak` = "(b) Peak", `End of Season` = "(c) End of Season", `Season Duration` = "(d) Season Duration"))) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.3, position = position_dodge(0.9)) +
  geom_text(aes(x = eco_state, y = emmean, label = .group), 
            vjust = -0.8, show.legend = FALSE, position = position_dodge(0.9)) +
  coord_flip() +
  ylab("Day of Year (Duration: Number of Days)") +
  xlab("") +
  scale_x_discrete(labels=c("Gravelly\nShrubland","Sandy\nShrub-invaded\nGrassland","Sandy\nShrubland")) +
  labs(col = "Rainfall") +
  scale_color_discrete(breaks=c('wet','average','dry')) +
  theme(panel.spacing.y = unit(0, "lines"))

scale_individual_facet_y_axes = function(plot, ylims) {
  init_scales_orig = plot$facet$init_scales
  
  init_scales_new = function(...) {
    r = init_scales_orig(...)
    # Extract the Y Scale Limits
    y = r$y
    # If this is not the y axis, then return the original values
    if(is.null(y)) return(r)
    # If these are the y axis limits, then we iterate over them, replacing them as specified by our ylims parameter
    for (i in seq(1, length(y))) {
      ylim = ylims[[i]]
      if(!is.null(ylim)) {
        y[[i]]$limits = ylim
      }
    }
    # Now we reattach the modified Y axis limit list to the original return object
    r$y = y
    return(r)
  }
  
  plot$facet$init_scales = init_scales_new
  
  return(plot)
}

ylims = list(c(25,225), c(100,300), c(200,400), c(75,275), c(25,225), c(100,300), c(200,400), c(75,275))
scale_individual_facet_y_axes(plot, ylims = ylims)
ggsave("figures/LSMeans 4-way interaction.jpg", height = 8, width = 10, units = "in", dpi = 300)

##########################################
# FIGURE: HLS PEAK EVI LEAST SQUARE MEANS
##########################################

hls_eco_state_rainfall_lsm <- hls.eco_state.rainfall.lsm %>%
  mutate(response = "Peak EVI") %>%
  mutate(.group = trimws(.group)) %>%
  mutate(rainfall_factor = factor(rainfall, levels = c("dry", "average", "wet")))

ggplot(hls_eco_state_rainfall_lsm, aes(x = eco_state, y = emmean, col = rainfall_factor)) +
  theme_bw() +
  geom_point(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), width=.3, position = position_dodge(0.9)) +
  geom_text(aes(x = eco_state, y = emmean, label = .group), 
            vjust = -0.8, show.legend = FALSE, position = position_dodge(0.9)) +
  coord_flip() +
  ylab("Peak EVI") +
  xlab("") +
  labs(col = "Rainfall") +
  scale_color_discrete(breaks=c('wet','average','dry'))

ggsave("figures/Peak EVI LSMeans EcoState-Rainfall interaction.jpg", height = 6, width = 8, units = "in", dpi = 300)