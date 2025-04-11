# Alberto Rovellini
# 03/31/2025
# This code takes ROMS-NPZ output for historical and projection runs from our GOACLIM NEP10K monthly indices by NMFS area and transforms them in input forcing files for Atlantis
# The goal is to produce a scalar from present days to apply in the projections

# there are many ways to do this:
# Force plankton time series (not doing that as it takes recalibration)
# Force mum scalars based on monthly and box-specific changes compared to historical period
# Get total trend over the time series and impose that as a scalar

# The most "thorough" approach would be forcing plankton, but that's a phylosophy choice in a way - and it would also result into a new model
# So passing that, I think that keeping things simple makes sense

# The most straightforward approach is to compute a change from present days to end of century and then apply it linearly
# Do it on production for all the plankton groups
# Do it equally across space as to simplify
# Many choices to make. What time frame? What time of the year? What water layer? Use the hindcast or not?
# Layer: we are talking about plankton production so what happens at the surface seems more relevant
# Time of the year: given the bloom dynamics in the GOA, spring seems sensible (though there is a lot of variability between groups)
# Time period: 10 years
# Hindcast: I'd say no. The hindcast has data, but the projection does not. So, if the goal is the relative change, raw historical / projection sounds best
# Start: 2015:2024
# End: 2090:2099

library(tidyverse)
library(lubridate)
library(zoo)
library(forecast)

# read data, use average instead of sum
nep_hist <- read.csv("data/NEP_10k_revised_indices/nep_sum_wb_hist_1000.csv")
nep_hist$simulation = "historical"

nep_ssp126 <- read.csv("data/NEP_10k_revised_indices/nep_sum_wb_ssp126_1000.csv")
nep_ssp126$simulation = "ssp126"

nep_ssp245 <- read.csv("data/NEP_10k_revised_indices/nep_sum_wb_ssp245_1000.csv")
nep_ssp245$simulation = "ssp245"

nep_ssp585 <- read.csv("data/NEP_10k_revised_indices/nep_sum_wb_ssp585_1000.csv")
nep_ssp585$simulation = "ssp585"

roms_avg_data <- do.call(rbind, list(nep_ssp126, nep_ssp245, nep_ssp585))

# what are the available names?
plankton_vars <- c("prod_Cop", "prod_Eup", "prod_MZL", "prod_MZS", "prod_PhL", "prod_PhS")

roms_avg_data <- roms_avg_data %>% 
  filter(varname %in% plankton_vars) %>%
  mutate(
    date = lubridate::as_date(date),
    month = lubridate::month(date),
    year = lubridate::year(date))

# view in time
roms_avg_data %>%
  filter(NMFS_AREA == "All", month == 7) %>%
  ggplot(aes(x = date, y = value, color = simulation)) +
  geom_line()+
  facet_wrap(depthclass~varname, scales = "free")

# mean of first 10 years
# mean of last 10 years
# do spring surface eup, diatom, cop
# get clim
roms_summary <- roms_avg_data %>%
  filter(NMFS_AREA == "All", # all areas on the shelf
         depthclass == "Surface", # surface only
         month %in% c(4:6)) %>% # spring quarter
  group_by(varname, simulation, year) %>%
  summarise(meanvar = mean(value))

# present day hindcast
roms_now <- roms_summary %>%
  filter(between(year, 2015,2024)) %>%
  group_by(varname, simulation) %>%
  summarise(clim = mean(meanvar)) %>%
  rename(now = clim)

# end of century
roms_eoc <- roms_summary %>%
  filter(between(year, 2090,2099)) %>%
  group_by(varname, simulation) %>%
  summarise(clim = mean(meanvar)) %>%
  rename(eoc = clim)

# combine
roms_change <- roms_now %>%
  left_join(roms_eoc) %>%
  mutate(change = eoc / now)

# view
roms_change %>%
  ggplot(aes(x = varname, y = change, fill = simulation))+
  geom_bar(stat = "identity", position = position_dodge())

# print
roms_change

# save
saveRDS(roms_change, "output/roms_change.RDS")

# #######################
# # Approach with time series decomposition
# # Not using this but keeping it for future reference
# # Load required libraries
# library(stats)
# library(ggplot2)
# library(dplyr)
# library(lubridate)
# library(trend) # For Mann-Kendall test
# 
# # Step 1: Data preparation
# ts_data <- roms_avg_data %>%
#   filter(varname == "prod_PhS",
#          NMFS_AREA == "All",
#          depthclass == "Surface",
#          simulation == "ssp585") %>%
#   # Ensure data is ordered chronologically
#   arrange(date) %>%
#   # Select relevant columns
#   select(date, value, year, month)
# 
# # Step 2: Convert to ts object
# # First check if you have complete monthly data (important for ts object creation)
# date_range <- seq(min(ts_data$date), max(ts_data$date), by = "month")
# if(length(date_range) != nrow(ts_data)) {
#   print("Warning: You have missing months in your data. Consider imputation.")
#   # Simple imputation example (you may want a more sophisticated approach)
#   all_dates <- data.frame(date = date_range)
#   ts_data <- full_join(all_dates, ts_data, by = "date") %>%
#     arrange(date) %>%
#     mutate(
#       year = year(date),
#       month = month(date)
#     )
#   # Simple linear interpolation for missing values
#   ts_data$value <- na.approx(ts_data$value, na.rm = FALSE)
# }
# 
# # Create ts object
# start_year <- min(ts_data$year)
# start_month <- min(ts_data$month[ts_data$year == start_year])
# ts_obj <- ts(ts_data$value, 
#              start = c(start_year, start_month), 
#              frequency = 12) # 12 for monthly data
# 
# # Step 3: Apply STL decomposition
# stl_result <- stl(ts_obj, s.window = "periodic", robust = TRUE)
# 
# # Step 4: Extract components
# trend <- stl_result$time.series[, "trend"]
# seasonal <- stl_result$time.series[, "seasonal"]
# remainder <- stl_result$time.series[, "remainder"]
# 
# # Step 5: Visualize decomposition
# plot(stl_result, main = "STL Decomposition")
# 
# # Step 6: Calculate trend change using 10-year windows
# # Modified to use 10 years (120 months) at each end
# first_period <- 1:120  # First 10 years of monthly data
# last_period <- (length(trend) - 119):length(trend)  # Last 10 years
# 
# first_decade_avg <- mean(trend[first_period], na.rm = TRUE)
# last_decade_avg <- mean(trend[last_period], na.rm = TRUE)
# 
# absolute_change <- last_decade_avg - first_decade_avg
# percent_change <- (absolute_change / first_decade_avg) * 100
# 
# # Step 7: Statistical significance testing
# mk_test <- mk.test(trend)
# sen_slope <- sens.slope(trend)
# 
# # Calculate annualized change (divide by number of years)
# years_span <- max(ts_data$year) - min(ts_data$year) + 1
# annual_change_rate <- absolute_change / years_span
# annual_percent_rate <- percent_change / years_span
# 
# # Print results
# cat("Trend Analysis Results (using 10-year averages):\n")
# cat("First decade average (", start_year, "-", start_year + 9, "): ", first_decade_avg, "\n")
# cat("Last decade average (", max(ts_data$year) - 9, "-", max(ts_data$year), "): ", last_decade_avg, "\n")
# cat("Total absolute change: ", absolute_change, " ", unique(ts_data$unit), "\n")
# cat("Total percent change: ", percent_change, "%\n")
# cat("Annual rate of change: ", annual_change_rate, " ", unique(ts_data$unit), " per year\n")
# cat("Annual percent rate: ", annual_percent_rate, "% per year\n")
# cat("Mann-Kendall test p-value: ", mk_test$p.value, "\n")
# cat("Sen's slope (median change per month): ", sen_slope$estimates, "\n")
# cat("Sen's slope (annual rate): ", sen_slope$estimates * 12, " per year\n")
# 
# # Optional: Add a visualization of the trend with 10-year averages highlighted
# years <- seq(start_year, max(ts_data$year), length.out = length(trend))
# trend_df <- data.frame(Year = years, Trend = trend)
# 
# ggplot(trend_df, aes(x = Year, y = Trend)) +
#   geom_line() +
#   geom_rect(aes(xmin = min(years), xmax = min(years) + 10, 
#                 ymin = -Inf, ymax = Inf), 
#             fill = "lightblue", alpha = 0.3) +
#   geom_rect(aes(xmin = max(years) - 10, xmax = max(years), 
#                 ymin = -Inf, ymax = Inf), 
#             fill = "lightblue", alpha = 0.3) +
#   geom_hline(yintercept = first_decade_avg, linetype = "dashed", color = "blue") +
#   geom_hline(yintercept = last_decade_avg, linetype = "dashed", color = "red") +
#   labs(title = "Long-term Trend with 10-year Average Windows Highlighted",
#        x = "Year", 
#        y = "Trend Component") +
#   theme_minimal()
