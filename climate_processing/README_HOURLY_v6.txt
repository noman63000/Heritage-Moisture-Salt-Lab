Mohenjo-daro Climate Forcing Package v6
Generated 28 September 2026

PURPOSE
This package contains the climate boundary forcing prepared for the secondary-data-driven Mohenjo-daro heat-moisture-salt simulation.

HISTORICAL / REFERENCE FORCING
- ERA5-Land 2010-2023, 122,712 hourly records.
- Model-ready file: Mohenjo_daro_ERA5Land_2010_2023_hourly_model_ready.csv.gz
- Grid point: 27.3 N, 68.1 E.
- Tiny negative de-accumulation artifacts in precipitation and shortwave radiation were clamped to zero.

FUTURE FORCING
- NASA NEX-GDDP-CMIP6 v2 daily forcing.
- Models: ACCESS-CM2, CanESM5, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0.
- Scenarios: SSP2-4.5 and SSP5-8.5.
- Periods: 2040-2060 and 2080-2100.
- 20 independent hourly forcing files, 3,681,312 total hourly rows.
- NEX grid point: 27.375 N, 68.125 E.

TEMPORAL DISAGGREGATION
Daily NEX values were converted to hourly forcing using same-month ERA5-Land analog-day profiles selected in quantile space. NEX daily mean/min/max temperature, mean RH, precipitation total, mean shortwave, mean longwave and mean wind speed are preserved to numerical precision.

Physical controls:
- Negative de-accumulation artifacts in precipitation/shortwave were set to zero and each affected day renormalized to the original NEX daily target.
- Rare scaled shortwave values above 1200 W/m2 were capped and the removed energy redistributed over daylight hours while preserving the daily mean. 440 of 153,388 future model-days reached this ceiling in the final forcing.
- Final QA: no negative physical flux values and no RH values outside 0-100%.

IMPORTANT LIMITATION
This is statistical temporal disaggregation. It does not independently project future sub-daily storm sequencing or hourly extremes. Hourly structure is inherited from historical ERA5-Land analogs. This limitation should be disclosed in the manuscript and tested in sensitivity analysis.

SIMULATION STRATEGY
Do not average the five GCM climate forcings before the nonlinear deterioration model. Run the historical/reference case and all 20 future cases independently, then summarize deterioration outputs across GCMs (for example median and IQR).

FILES
- Mohenjo_daro_Simulation_Input_Data_v6_HOURLY_READY.xlsx: master input/provenance workbook.
- hourly_future_v6/: 20 future hourly forcing CSV.gz files plus QA, checksums, analog mapping and metadata.
- Mohenjo_daro_ERA5Land_2010_2023_hourly_model_ready.csv.gz: historical/reference hourly forcing.
- Mohenjo_daro_NEX_CMIP6_daily_all_models.csv.gz: canonical daily future forcing used as the exact daily target.
- Mohenjo_daro_ERA5Land_2010_2023_raw_merged.csv.gz: raw historical source used for temporal profiles.
- processing_scripts/: scripts used to generate, clean and validate the hourly forcing.

NEXT MODEL BLOCKER
Climate forcing is ready. The numerical heat-moisture-salt model should not be run as a publishable study until the remaining material/transport parameters are locked from defensible sources: thermal conductivity, specific heat, vapour transport, liquid/capillary transport, sorption/moisture-storage behavior, salt diffusion/dispersion and the exact crystallization formulation/kinetic parameters (or a justified equilibrium alternative).
