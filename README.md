# SeaLevelExtrapolation_2024
This repository contains the codes and input data to reproduce the results from the following paper:

A.S. Bellas-Manley, R.S. Nerem, & B.D. Hamlington: Extrapolation of the Satellite Altimeter Record to Understand Regional Variations in Future Sea Level Change, JGR:Oceans, submitted

## User Manual 

The scripts in this directory can be used to perform a quadratic extrapolation of sea level based on a gridded data set like the MEaSUREs Sea Surface Height data product. The tasks performed include (i) computing the regional timeseries by averaging the gridded data in regions specified by a 5-point polygon, (ii) computing the rate and acceleration of the regional timeseries data, (iii) computing the errors associated with serially correlated formal errors, GIA, and measurement errors, (iv) extrapolating the preferred pair of rate & acceleration as well as the confidence limits based on the errors from (iii), and finally (v) plotting the results in comparison with the AR6 projections.

Software required to run the scripts: a Python installation with the following packages installed: `SciPy`, `NumPy`, `matplotlib`, `netCDF`.

### Data 

Data required to run the scripts can be downloaded from the sources below, or is included in the directory data/

MEaSUREs Gridded Sea Surface Height Anomalies: https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL2205

AR6 medium confidence sea level projections: https://zenodo.org/records/6382554

### Steps
1. Run Assemble_ssh_timeseries_file.py
    - This script reads the individual MEaSUREs data files and writes to a single file after downsmapling spatially and averaging monthly
1.	Run DefineRegions.py
    - This script contains the coordinates that define a five-point polygon containing the desired region
    - For a lower order polygon region, just repeat one or two of the coordinates
    - First coordinate must be the lower left corner and list coordinates counter-clockwise
1.	Run IsolateRegionalTimeseries_ComputeRateAccel.py
    - Calls a function IsolateTimeseries_5PointPolygonRegions to 
        - Compute the regional timeseries by averaging the MEaSUREs data within defined regions
        - Save results to file 
    - Calls a function ComputeRateAccel_RegionalTimeseries to 
        - compute the rate, acceleration, annual, semi-annual, and formal errors of the regional timeseries, 
        -save results to file
    - Plots the results 
1.	Run ExtrapolateRegionalTimeseries_PlotwithAR6multiSSP.py
    - Calls a function ConstructRateAccelEnsemble to 
        - Compute the regional error associated with GIA, serially correlated formal errors, and measurement errors
        - construct an ensemble of pairs of rate and acceleration based on the errors
        - return the distribution of rates and accelerations
    - Calls a function ComputeExtrapolations_Preferred_and_ConfidenceLimits to 
        - extrapolate the preferred regional trajectory of sea level, and the confidence limits,
    - Calls a function IsolateRegionalTimeseriesAR6_5PointPolygon to
        - Compute the regional timeseries by averaging the AR6 medium confidence total sea level projections within defined regions
    - Loads the regional AR6 timeseries data, 
    - Calls a function SetExtrapolationReferenceYear to 
    - Reference all timeseries to a specified year, and
    - Plots the results
