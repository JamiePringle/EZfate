# EZfate -- global estimates of connectivity in and near the coastal ocean. 

## Software in R for obtaining and analyzing precomputed Lagrangian Pathways computed with currents from the Mercator 1/12th degree ocean model

This repository contains code for analyzing precomputed Lagrangian pathways in the coastal ocean from the Mercator System For Global Ocean Physical Analysis At 1/12° as described [in this data sheet](https://www.mercator-ocean.eu/wp-content/uploads/2017/02/SYSTEM-sheet-_PSY4V3R1_2017.pdf). The particle tracking was made with oceanParcels, a python package described in at [oceanparcels.org](https://oceanparcels.org/).

For instructions on how to download and use this code, go to [the instructions in the docs directory](https://jamiepringle.github.io/EZfate/)

By September 2023, the data and software will produce estimates of connectivity for indvidual months and climatological months for particles released at 1m, 10m, 20m, and 40m depths and drifting for 2, 4, 8 ... 60 days for all locations within 100km of the coast or depths of less than 500m. The particles can either be fixed to their initial depth or allowed to advect to different depths by the grid-scale vertical velocity. In the future more depths will be added, along with the vertical dispersal of particles by turbulence as computed by the ocean model. 

Monthly and climatological data are formed from twice-a-day particle releases at all locations. 

Pathways have been calculated for the years 2007-2020 for all of the depths above for the coasts of North and South America. 

Global pathway calculation is underway now for years 2007-2022. The status of depth of release and depth keeping behavior is given below. Blank indicates that it has not yet been done. "Fixed Depth Behavior" indicates particles are only advected in horizontal directions by ocean currents, "Variable Depth Behavior" indicates that they are advected by both horizontal and vertical behavior. 

|Depth|Fixed Depth Behavior|Variable Depth Behavior|
|-----|--------------------|-----------------------|
| 1m |  | Done|
|10m | Done| |
|20m | | |
|40m | | |

This work is supported by NSF project OCE 1947954.  

# Please direct questions to James Pringle, University of New Hampshire, jpringle@unh.edu. 
