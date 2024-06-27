# SeaLevelExtrapolation_2024
This repository contains the codes and input data to reproduce the results from the following paper:

A.S. Bellas-Manley, R.S. Nerem, & B.D. Hamlington: Extrapolation of the Satellite Altimeter Record to Understand Regional Variations in Future Sea Level Change, JGR:Oceans, submitted

## User Manual 

The scripts in this directory can be used to perform a quadratic extrapolation of sea level based on a gridded data set like the MEaSUREs Sea Surface Height data product. The tasks performed include (i) computing the regional timeseries by averaging the gridded data in regions specified by a 5-point polygon, (ii) computing the rate and acceleration of the regional timeseries data, (iii) computing the errors associated with serially correlated formal errors, GIA, and measurement errors, (iv) extrapolating the preferred pair of rate & acceleration as well as the confidence limits based on the errors from (iii), and finally (v) plotting the results in comparison with the AR6 projections.

Software required to run the scripts: a Python installation with the following packages installed: `SciPy`, `NumPy`, `matplotlib`, `netCDF`.

All data required to run the scripts is included in the directory data/

### Steps
1.	Run DefineRegions.py
    -This script contains the coordinates that define a five-point polygon containing the desired region
    -For a lower order polygon region, just repeat one or two of the coordinates
    -First coordinate must be the lower left corner and list coordinates counter-clockwise
1.	Run IsolateRegionalTimeseries_ComputeRateAccel.py
  1.	Calls a function IsolateTimeseries_5PointPolygonRegions to 
     1.	Compute the regional timeseries by averaging the MEaSUREs data within defined regions
     1.	Save results to file 
  b.	Calls a function ComputeRateAccel_RegionalTimeseries to 
    i.	compute the rate, acceleration, annual, semi-annual, and formal errors of the regional timeseries, 
    ii.	save results to file
  c.	Plots the results 
1.	Run ExtrapolateRegionalTimeseries_PlotwithAR6multiSSP.py
  a.	Calls a function ConstructRateAccelEnsemble to 
    i.	Compute the regional error associated with GIA, serially correlated formal errors, and measurement errors
    ii.	construct an ensemble of pairs of rate and acceleration based on the errors
    iii.	return the distribution of rates and accelerations
  b.	Calls a function ComputeExtrapolations_Preferred_and_ConfidenceLimits to 
    i.	extrapolate the preferred regional trajectory of sea level, and the confidence limits,
  c.	Calls a function IsolateRegionalTimeseriesAR6_5PointPolygon to
    i.	Compute the regional timeseries by averaging the AR6 medium confidence total sea level projections within defined regions
  d.	Loads the regional AR6 timeseries data, 
  e.	Calls a function SetExtrapolationReferenceYear to 
  i.	Reference all timeseries to a specified year, and
  f.	Plots the results
