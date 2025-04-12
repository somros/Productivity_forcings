# Alberto Rovellini
# 4/12/2025
# Code to write external scaling input for plankton growth rates
# ratios are created in plankto_scalars.R

# See https://confluence.csiro.au/display/Atlantis/External+Mortality%2C+Growth+and+Recruitment+Scaling

# this exercise is hinged on the following mapping from NPZ to Atlantis GOA:
# time-averaged small coastal copepod concentration <- Mesozooplankton
# time-averaged euphausiid concentration <- Euphausiids
# time-averaged large microzooplankton concentration <- Microzooplankton
# time-averaged small microzooplankton concentration <- Microzooplankton
# time-averaged neocalanus spp. concentration <- Mesozooplankton
# time-averaged large phytoplankton concentration <- Diatoms
# time-averaged small phytoplankton concentration <- Picophytoplankton

library(tidync)
library(dplyr)
library(tidyr)
library(ncdf4)
library(RNetCDF)
library(lubridate)
library(ggplot2)

options(digits=6)

# open ratios
ratios <- readRDS("output/roms_change_lm.RDS")

# Make nc scalar file --------------------------------------------
nbox <- 109
nlayer <- 7
t_units <- "seconds since 1990-01-01 00:00:00 -9" # for consistency with the init.nc files
seconds_timestep <- 43200 * 2 * 365
this_geometry <- "GOA_WGS84_V4_final.bgm"

# time array
# 30 years burn-in
# then start 2020 and do 2020-2100
burnin <- 30
simtime <- 80
nyr <- burnin + simtime
time_array <- seq(0, seconds_timestep * nyr, seconds_timestep)
ntime <- length(time_array)

ssp_vars <- unique(ratios$simulation)
code_vars <- unique(ratios$Code)

for(i in 1:length(ssp_vars)){
  ssp <- ssp_vars[i]
  
  nc_name <- paste0('output/scalar_proj_ROMS_', ssp, '.nc')
  nc_file <- create.nc(nc_name)
  
  # set and write dimensions and global attributes of the ncfile
  dim.def.nc(nc_file, "t", unlim=TRUE)
  dim.def.nc(nc_file, "b", nbox) # manual
  dim.def.nc(nc_file, "z", nlayer) # manual
  
  var.def.nc(nc_file, "t", "NC_DOUBLE", "t")
  
  att.put.nc(nc_file, "t", "units", "NC_CHAR", t_units)
  att.put.nc(nc_file, "t", "dt", "NC_DOUBLE", seconds_timestep)
  att.put.nc(nc_file, "NC_GLOBAL", "title", "NC_CHAR", "External scalar for plankton growth")
  att.put.nc(nc_file, "NC_GLOBAL", "geometry", "NC_CHAR", this_geometry)
  att.put.nc(nc_file, "NC_GLOBAL", "parameters", "NC_CHAR", " ")
  
  var.put.nc(nc_file, "t", time_array)
  
  # now handle the groups
  for(j in 1:length(code_vars)){
    var <- code_vars[j]
    
    target <- ratios %>% filter(simulation == ssp, Code == var) %>% pull(relative_change)
    ramp <- seq(from = 1, to = target, length.out = simtime)
    
    # Initialize the array with NAs
    scalar_array <- array(NA, dim = c(nlayer, nbox, nyr+1)) # you need the +1 because the file needs to start from t0
    
    # Fill first slices with 1s
    scalar_array[,,1:burnin+1] <- 1
    
    # Fill remaining slices with values from the vector
    for (k in 1:simtime) {
      scalar_array[,,burnin+1+k] <- ramp[k]
    }
    
    varlab <- paste0(var,"growth_ff")
    var.def.nc(nc_file, varlab, "NC_DOUBLE", c("z","b","t"))
    att.put.nc(nc_file, varlab, "_FillValue", "NC_DOUBLE", 1)
    var.put.nc(nc_file, varlab, scalar_array)
    
  }

  close.nc(nc_file)
  
}
