
### EZfate and python, a quick guide

The EZfate package on GitHub is written mostly in R, but most of the data was created in python, and all the anlysis can be performed in python. The data structures are available in zarr, which is easilly read by python. The following text describes the underlying data structures that EZfate provides, and documents some routines I have provided to get started using EZfate data in python. 

All of this description will make more sense if you spend some time reading at least the introduction to EZfate given at its GitHub page [here](https://github.com/JamiePringle/EZfate). 

#### Python packages you need

The codes I provide expect that the following python modules exist on your system. They are all easily installed with conda:

* numpy
* matplotlib
* zarr
* netCDF4

#### Data structures
EZfate provides two main data structures, the forward connectivity `E` and the backwards connectivity `Etranspose`. These are best thought of as sparce matrices that define where all the particles released at one location go (`E`) and where all particles released at one location came from (`Etranspose`). The two data structures are stored as Zarr files which are downloaded by `getEZfateFromOSN.py`, described below. They have been precomputed globally for all months from 2007 to last year, as is described on the EZfate GitHub page.  Climatological connectivity is also provided. 

All locations are given as indices into the Mercator Ocean Model GLORYS 1/12th degree grid. This was done to conserve space and bandwidth. Instructions for converting this into latitude and longitude is given below. 

`E` contains 6 columns of data: 

* `nxFrom` and `nyFrom`, integer vectors of the locations where particles are released from
* `nxTo, nyTo`, a vector of variable length vectors which gives the locations where particles released from (`nxFrom,nyFrom`) end up.
* `numTo`, a vector of variable length vectors which give the number of particles which end up in each (`nxTo,nyTo`)
* `numLaunched`, a vector of integers of the total number of drifters launched from each (nxFrom,nyFrom) point that are included in this data. This number can change a lot from point to point if we only include tracks that start and end in the habitat `E` – this is an option in the creation of `E`

`Etranspose` contains a similar 6 columns of data, and is nearly symmetric with `E` with "From" and "To" switched:

* `nxTo` and `nyTo`, integer vectors of the locations where arrive after some time
* `nxFrom, nyFrom`, a vector of variable length vectors which gives the locations where particles which arrived at (`nxTo,nyTo`) started from.
* `numFrom`, a vector of variable length vectors which give the number of particles which end up in each (`nxFrom,nyFrom`)
* `cumSumPropFrom`, which I will neglect to define now.

#### (nx,ny) to (longitude,latitude)

As part of the data distributed, the model grid for the GLORYS model is provided. It can be fetched with the provided `getEZfateFromOSN.py` module, and used to convert vectors of `nx,ny` to longitude and latitude with the following code. 
```python
import getEZfateFromOSN
import netCDF4 as nc

#get the file which contains the bathymetry and grid data for the mercator ocean model and for
#the connectivity files
maskFile='EZfateData/EZfateFiles/ext-PSY4V3R1_mesh_zgr.nc'
maskFile=getEZfateFromOSN.getFileFromOSN(maskFile)
with nc.Dataset(maskFile,'r') as gridData:
    lonMat=gridData['nav_lon'][:]
    latMat=gridData['nav_lat'][:]

#write a function that takes a vector of nx grid locations and ny
#locations and returns vectors of longitude and latitude.
def nxny2lonlat(nxVec,nyVec):
    lonVec=array([lonMat[p] for p in zip(nyVec,nxVec)])
    latVec=array([latMat[p] for p in zip(nyVec,nxVec)])
    return lonVec,latVec
```

#### Example files

Soon I will provide more documentation, but for now I will provide three helpful (I hope) example python scripts to show how to use this connectivity data. 

* `00_makeConnectivityMatrices_trimByDistance_climatologicalData.py` creates the connectivity matrices `E` and `Etranspose` for climatological months and for arbitrary regions (called 'habitats') for a given drift duration in days (currently, EZfate includes all even durations from 2 to 60 days)
* `00_makeConnectivityMatrices_trimByDistance_singleYear.py` creates the connectivity matrices `E` and `Etranspose` for months in a single year and for arbitrary regions (called 'habitats') for a given drift duration in days (currently, EZfate includes all even durations from 2 to 60 days)
* `01_plot_connectivty_and_transpose_entireHabitabitat.py` plots all starting (`From`) and ending (`To`) locations
* `02_plotConnectivity_to_and_from.py` finds a number of points in the habitat close to some arbitrary latitude and longitude, and plots both where all Lagrangian particles leaving from that point go and where all particles reaching that point came from

The graphics produced are crude because I did not want to add cartopy or equivalent as a dependency. 

#### Usefull python functions
These call several important utility modules in which much of the real work is done. In the future I will document the many useful functions in them in this document, but they are well documented in the files themselves. Below is a subset of the most useful functions they provide. 

* `getEZfateFromOSN.py` gets the EZfate files from the Open Storage Network s3 bucket
  * `getFileFromOSN(pathName)` gets a file from the bucket and saves it in a local directory
* `makeConnectivityModule.py` mostly is used to take Lagrangian data from oceanParcels and convert it into connectivity matrices, but has some useful utility functions I use.
  * `makeEmptyConnectivity()` makes a new `E` or `Etranspose` zarr file
  * `getGridDistanceDict()` gets the distance of each nx,ny point from land
  * `getGridDepthDict()` gets the depth of each nx,ny point
  * `trimConnectivity()` trims the `E` data to some habitat by depth or distance
  * `trimConnectivity_byPoly()` trims the `E` data to inside or outside some polygon
  * `combineConnectivity()` merges two `E` structures
* `popMatrixManipulate.py` contains many useful utility functions that manipulate connectivity matrices, like transposing, joining, and normalizing them.
  * `invertConMatrix()` converts a `E` to `Etranspose`

  