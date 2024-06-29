# EZfate -- global estimates of connectivity in and near the coastal ocean. 

## Software in R for obtaining and analyzing precomputed Lagrangian Pathways computed with currents from the Mercator 1/12th degree ocean model

<p align="center">
  <img src="https://jamiepringle.github.io/EZfate/twoStarts_EastFL_MxQR_AprMayJun.gif"><br>
  <em>Climatological dispersal of particles released from Quintana Roo, Mexico, and the East Coast of Florida, USA in April, May and June, and allowed to disperse for 18 days. Particles released at 1m, can vary in depth subsequently. Colors represent log of fraction of particles released which reach a point. Code to produce this plot is in code directory, and is in two files:</em> <tt>animate_DispersalData_from_two_Regions.Rmd</tt> <em>and</em> <tt>save_DispersalData_fromRegion.R</tt> in <tt>code/code_for_README.md_graphics/</tt>.
</p> 

This repository contains code for analyzing precomputed Lagrangian pathways in the coastal ocean from the Mercator System For Global Ocean Physical Analysis At 1/12° as described [in this data sheet](https://www.mercator-ocean.eu/wp-content/uploads/2017/02/SYSTEM-sheet-_PSY4V3R1_2017.pdf) and in [Lellouche et al (2018)](https://os.copernicus.org/articles/14/1093/2018/). The particle tracking was made with oceanParcels, a python package described in at [oceanparcels.org](https://oceanparcels.org/).

For instructions on how to download and use this code, go to [the instructions in the docs directory](https://jamiepringle.github.io/EZfate/) . If you wish to read the data into python, contact me and I can show you how to access the Zarr files with the connectivity data. Note well: if you have problems downloading the data, please update to the latest version of EZfate.

The data and software will produce estimates of connectivity for indvidual months and climatological months for particles released at 1m, 10m, 20m, and 40m depths and drifting for 2, 4, 8 ... 60 days for all locations within 100km of the coast or depths of less than 500m. The particles can either be fixed to their initial depth or allowed to advect to different depths by the grid-scale vertical velocity. In the future more depths will be added, along with the vertical dispersal of particles by turbulence as computed by the ocean model. 

Monthly data for each year and climatological data for all years run are formed from twice-a-day particle releases at all locations. *Note Well:* Each year at arround June, the climatological data will change as I add the prior years data. If you want a fixed set of years included in the climatological average, either save the data you download, or download the years and regions you desire and form your own climatology (both are described in the documentation).

Pathways have been calculated for the years 2007-2023 for all of the depths above for the coasts of North and South America. In the future, the prior year will be added by June of the following year. 

 The status of depth of release and depth keeping behavior is given below. "Fixed Depth Behavior" indicates particles are only advected in horizontal directions by ocean currents and remain at their starting depth. "Variable Depth Behavior" indicates that they are advected by both horizontal and vertical currents and can vary in depth. In the table, the last year which has been completed is given. 

|Depth|Fixed Depth Behavior|Variable Depth Behavior|
|-----|--------------------|-----------------------|
| 1m | 2023 | 2023 |
|10m | 2023 | 2023 |
|20m | 2023 | 2023 |
|40m | 2023 | 2023 |

This work is supported by NSF project OCE 1947954.  

**_Python software_**
If you wish to analyze the data using python, contact me, and I shall share my beta-version python code. The code is feature complete (it is used to calculate the R data shared above), but has not yet been fully documented. 

**_Suggested Citation_**
> Pringle, J.M. (2023). EZfate, a tool for estimating larval connectivity in the global coastal ocean [Computer software].
> [![DOI](https://zenodo.org/badge/569445832.svg)](https://zenodo.org/doi/10.5281/zenodo.10214924)

**_Suggested Acknowledgment_**

> A portion of this work used code provided by Jamie Pringle's EZfate package (https://zenodo.org/doi/10.5281/zenodo.10214924).

# Please direct questions to James Pringle, University of New Hampshire, jpringle@unh.edu. 
