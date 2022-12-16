# WARNING-- THIS CODE IS ACTIVELY BEING DEVELOPED. PLEASE DON'T USE FOR ANYTHING IMPORTANT, OR COMPLAIN IF ANYTHING YOU DO WITH IT BREAKS. Please direct questions to James Pringle, University of New Hampshire, jpringle@unh.edu


# EZfate -- global (soon) estimates of connectivity in and near the coastal ocean. 

## Software for analysis of precomputed Lagrangian Pathways from Mercator 1/12th degree ocean model

This repository contains code for analyzing precomputed Lagrangian pathways in the coastal ocean from the System For Global Ocean Physical Analysis At 1/12° as described [in this data sheet](https://www.mercator-ocean.eu/wp-content/uploads/2017/02/SYSTEM-sheet-_PSY4V3R1_2017.pdf). The particle tracking was made with oceanParcels, a python package described in at [oceanparcels.org](https://oceanparcels.org/).

For instructions on how to download and use this code, go to [the instructions in the docs directory](https://jamiepringle.github.io/EZfate/)

Currently, the data and software will produce estimates of connectivity for indvidual months and climatological months for particles released at 1m, 20m, and 40m depths and drifting for 2, 4, 8 ... 60 days. The particles can either be fixed to their initial depth or allowed to advect to different depths. In the future more depths will be added, along with the vertical dispersal of particles by turbulence as computed by the ocean model. 

The precomputed pathways were made with software from the OceanParcels project (https://oceanparcels.org).  It currently includes the years 2007-2020, but will extended to include through 2022 by the end of the first quarter of 2023.  This work is supported by NSF project OCE 1947954.  

The current available data include the coastal oceans of North and South America for all 1/12° separated points from the coast to either depths up to 500m deep or 100km from the coast, whichever is greater. Future data releases will be global. 

