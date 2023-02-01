# Please direct questions to James Pringle, University of New Hampshire, jpringle@unh.edu. This is beta software for now, so anticipate evolution of the code and product. A first official release is planned for early summer 2023.


# EZfate -- global (soon) estimates of connectivity in and near the coastal ocean. 

## Software in R for obtaining and analyzing precomputed Lagrangian Pathways computed with currents from the Mercator 1/12th degree ocean model

This repository contains code for analyzing precomputed Lagrangian pathways in the coastal ocean from the System For Global Ocean Physical Analysis At 1/12° as described [in this data sheet](https://www.mercator-ocean.eu/wp-content/uploads/2017/02/SYSTEM-sheet-_PSY4V3R1_2017.pdf). The particle tracking was made with oceanParcels, a python package described in at [oceanparcels.org](https://oceanparcels.org/).

For instructions on how to download and use this code, go to [the instructions in the docs directory](https://jamiepringle.github.io/EZfate/)

Currently, the data and software will produce estimates of connectivity for indvidual months and climatological months for particles released at 1m, 10m, 20m, and 40m depths and drifting for 2, 4, 8 ... 60 days. The particles can either be fixed to their initial depth or allowed to advect to different depths by the grid-scale vertical velocity. In the future more depths will be added, along with the vertical dispersal of particles by turbulence as computed by the ocean model. 

Pathways have been calculated for the years 2007-2020, but will extended to include through 2022 by the end of the first quarter of 2023.  

The current available data include the coastal oceans of North and South America for all 1/12° separated points from the coast to either depths up to 500m deep or 100km from the coast, whichever is greater. Future data releases will be global (Summer 2023). 

This work is supported by NSF project OCE 1947954.  

