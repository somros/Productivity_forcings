# Plankton Growth Scaling Factors for Atlantis GOA

This repository contains R scripts to generate external scaling inputs for plankton growth rates in the Atlantis ecosystem model for the Gulf of Alaska.

## Overview

These scripts calculate projected changes to plankton production based on ROMS-NPZ (Regional Ocean Modeling System - Nutrient-Phytoplankton-Zooplankton) output from projection climate runs. They create NetCDF scalar files that can be used as forcing in Atlantis GOA model projections.

## Contents

- `plankton_scalars.R`: Calculates relative changes in plankton production from ROMS-NPZ output
- `make_scalar.R`: Creates NetCDF scalar files for use in Atlantis model projections

## Workflow

1. `plankton_scalars.R` processes ROMS-NPZ output data to:
   - Extract and aggregate plankton production variables by Atlantis functional groups
   - Map NPZ plankton groups to Atlantis functional groups
   - Calculate relative changes in production using linear model fits from 2015 to 2100
   - Generate projections for different climate scenarios (SSP126, SSP245, SSP585)
   - Save results to `output/roms_change_lm.RDS`

2. `make_scalar.R` uses the calculated ratios to:
   - Create NetCDF scalar files for each climate scenario
   - Apply scaling factors to plankton growth rates
   - Generate linear ramps of change over the simulation period

## Data Mapping

The scripts use the following mapping from NPZ to Atlantis GOA:
- Small coastal copepod → Mesozooplankton (ZM)
- Euphausiid → Euphausiids (EUP)
- Large microzooplankton → Microzooplankton (ZS)
- Small microzooplankton → Microzooplankton (ZS)
- Neocalanus spp. → Mesozooplankton (ZM)
- Large phytoplankton → Diatoms (PL)
- Small phytoplankton → Picophytoplankton (PS)

## Usage

1. Ensure all required R packages are installed:
   ```R
   install.packages(c("tidync", "dplyr", "tidyr", "ncdf4", "RNetCDF", 
                     "lubridate", "ggplot2", "tidyverse", "zoo", "forecast"))
   ```

2. Run `plankton_scalars.R` first to generate the relative change values

3. Run `make_scalar.R` to create the NetCDF scalar files for Atlantis

## Output

The scripts generate:
- `output/roms_change_lm.RDS`: R data object containing relative changes in plankton production
- NetCDF scalar files: `output/scalar_proj_ROMS_[ssp]_NoBurnin.nc` for each climate scenario (SSP126, SSP245, SSP585)

## References

For more information on external scaling in Atlantis, see:
- [Atlantis External Mortality, Growth and Recruitment Scaling](https://confluence.csiro.au/display/Atlantis/External+Mortality%2C+Growth+and+Recruitment+Scaling)

## Notes

These scalars are needed because simply forcing increased temperature and salinity from ROMS forcings does not lead to realistic responses. We are missing the half of the story relative to productivity, which we have found is really important for this region and model (Rovellini et al. 2024).
One option in this sense would be forcing the nutrients and plankton in Atlantis with the ROMS-NPZ output, which some models do (e.g. NEUS, AMPS, etc.). This would require extensive recalibration (other than being more of a philosophical choice), and is not attainable at this time.
So, one option easier to plug in is to apply external forcings to the `mum` parameter of the relevant plankton groups. There are lots of ways to do this, but a very simple one to start is the approach outlined here. The high spatio-temporal aggregation is deliberate to keep things simple and acknowledge that these are more scenarios than realistic future conditions of productivity. 
Overall changes in production variables from NPZ models are somewhat modest, with the biggest change being a -25% for Euphausiids under ssp585.

## Author

Alberto Rovellini (2025)
arovel@uw.edu
