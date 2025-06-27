from pylab import *
from numpy import *
import zarr
from collections import defaultdict,Counter
import os
import shutil
import getEZfateFromOSN
import netCDF4 as nc

#This code shows how to plot all the starting points in a habitat, and where all of the particles go. 

#first, define the names of the connectivty and transpose files you wish to plot
if True:
    #for plotting climatological data
    baseName='CmaenasHab_depth1_minPLD40_maxPLD40_months5_to_6.zarr'
else:
    #plot a single year
    baseName='CmaenasHab_year2020_depth1_minPLD40_maxPLD40_months5_to_5.zarr'
    
#and then get the files
EfileName='connectivity/E_'+baseName
EtransposeFileName='connectivity/Etranspose_'+baseName
E=zarr.open(EfileName,'r')
Etranspose=zarr.open(EtransposeFileName,'r')

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

#now in figure 1, plot habitat locations in lon,lat space and plot all
#places where particles from those locations end up
figure(1);clf();style.use('ggplot')

#first, find all points that larvae can leave from
lonFrom,latFrom=nxny2lonlat(E.nxFrom,E.nyFrom)
plot(lonFrom,latFrom,'k*',markersize=5,zorder=10)

#now, for each point larvae can leave from, there are many points they
#can go to. So for each lonFrom and latFrom, there are many lonTo and
#latTo points. So while E.nxFrom[3] returns a single point, E.nxTo[3]
#and E.nyTo[3] will return a list of nx and ny coordinantes, and
#E.numTo[3] will return an integer giving how many drifters will go to
#each point. So a default dictionary will be used to count up how many
#particles go to each point. This double loop is certainly not the most
#efficient way to do this. But to make it a little faster, lets load the
#To variables into memory
EnxTo=E.nxTo[:]; EnyTo=E.nyTo[:]; EnumTo=E.numTo[:]
totalTo=defaultdict(lambda :0)
for nPoint in range(len(E.nxFrom)):
    for pNum in zip(EnxTo[nPoint],EnyTo[nPoint],EnumTo[nPoint]):
        p=pNum[:2]
        totalTo[p]=totalTo[p]+pNum[2]

nxToVecTotal=[p[0] for p in totalTo]
nyToVecTotal=[p[1] for p in totalTo]
numTotal=[totalTo[p] for p in totalTo]

#now make a scatter plot of where they end up
lonToTotal,latToTotal=nxny2lonlat(nxToVecTotal,nyToVecTotal)
scatter(lonToTotal,latToTotal,s=5,c=numTotal,vmin=0,vmax=400,zorder=5)
colorbar()

title('habitat')
draw()
show()
