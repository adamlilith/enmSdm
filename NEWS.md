enmSdm 0.5.2.7 2020-10-05
Added new methods and a fail-safe to interpolateRasters()

enmSdm 0.5.2.6 2020-10-03
Updated bioticVelocity() to use rasters, rather than arrays

enmSdm 0.5.2.4 2020-09-29
Fixed bug in interpolateRasters() when only one output raster is being interpolated

enmSdm 0.5.2.3 2020-09-24
Cell-weighted calculation of quantile lat/long in bioticVelocity()

enmSdm 0.5.2.3 2020-09-23
Remove ability bioticVelocity() to use "pophist" objects
New, straightforward examples in bioticVelocity()
Fixed "quant" and "Quants" velocities in bioticVelocity()

enmSdm 0.5.2.2
Graceful catch of non-converged/insufficient models in trainBrt()
Sensical names in output of bioticVelocity()

enmSdm 0.5.2.1
Consistent output across models for summaryByCrossValid()
Clarified help for bioticVelocity()

enmSdm 0.5.2.0
trainNs() can now return a list of all models
trainNs() now cycles through degrees of freedom
Can now use natural splines (trainNS()) with cross validation functions
Fixed issues with trainByCrossValid()

enmSdm 0.5.1.9
enmSdm 0.5.1.8
Create rastWithSquareCells()
Fix squareRasterCells(); renamed to squareRastCells()
Remove export of localSpatialCorrForValues() because of instability issue

enmSdm 0.5.1.7
Add squareRasterCells()

enmSdm 0.5.1.6
coordPrecision() handles coordinates in DMS format

enmSdm 0.5.1.5
Add Robinson projection to getCRS()

enmSdm 0.5.1.3
Recompiled for R 4.0.0

enmSdm 0.5.1.2
trainBrt() parallelized
trainMaxEnt() parallelized

enmSdm 0.5.0.4
trainByCrossValid() handles non-converged models

enmSdm 0.5.0.3
bioticVelocity() can correct for non-shared land mass
Fix bug in getCRS()

enmSdm 0.5.0.2
Fix bug in bioticVelocity()

enmSdm 0.5.0.1
Fix bug in spatialCorrForValues()

enmSdm 0.5.0.0
Fix bug in predictMaxEnt()

enmSdm 0.4.0.5
trainMaxEnt() and trainMaxNet() internals harmonized
Made sensible examples for trainXYZ() functions

enmSdm 0.4.0.3
Added function makeCRS() for making custom coordinate reference systems

enmSdm 0.4.0.2
Added function decimalToDms() for converting decimal coordinates to degrees-minutes-seconds format
Fixed bug in spatialCorrForPoints() when supplying distance bins in km

enmSdm 0.4.0.1
Fixed bug in coordPrecision() when supplied value was NA.

enmSdm 0.4.0.0
Fixed bug in elimCellDups() when there were no cell duplicates.