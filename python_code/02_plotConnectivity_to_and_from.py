from pylab import *
from numpy import *
import zarr
from collections import defaultdict,Counter
import os
import shutil
import getEZfateFromOSN
import netCDF4 as nc

#this code finds Npnts release points closest to a latitude/longitude
#point, and illustrates where particles released from those points go,
#and where particles that return to those points came from.

Npnts=20 #number of points to examine
lonP=-53.11; latP=48.685 #point to examine

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
#the transpose files
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

#now find the Npoint release points in the habitat closest to lonP,latP
#it would be more efficient and accurate to use a plane-sailing formula,
#but this should be good enough, is simple, and requires no external libraries
lonFrom,latFrom=nxny2lonlat(E.nxFrom,E.nyFrom) #lon,lat of release (From) points
scaleDistSq=((lonFrom-lonP)*cos(deg2rad(latFrom)))**2+(latFrom-latP)**2 #in weird units
closePoints=argsort(scaleDistSq)[:Npnts]


#now in figure 1, plot where points released from the Npoint points
#nearest to lonP,latP end up
figure(1);clf();style.use('ggplot')

#plot the location of the habitat, here the location of all lonFrom,latFrom points
plot(lonFrom,latFrom,'k,',markersize=5,zorder=10)
plot(lonFrom[closePoints],latFrom[closePoints],'r*',zorder=11)

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
for nPoint in closePoints:
    for pNum in zip(EnxTo[nPoint],EnyTo[nPoint],EnumTo[nPoint]):
        p=pNum[:2]
        totalTo[p]=totalTo[p]+pNum[2]

nxToVecTotal=[p[0] for p in totalTo]
nyToVecTotal=[p[1] for p in totalTo]
numTotal=[totalTo[p] for p in totalTo]

#now make a scatter plot of where they end up
lonToTotal,latToTotal=nxny2lonlat(nxToVecTotal,nyToVecTotal)
scatter(lonToTotal,latToTotal,s=5,c=numTotal,vmin=0,vmax=100,zorder=12)
colorbar()

title('where particles released from red stars ended up\n'+
      'colors represent count of particles')
draw()
show()

#====================================================================
#In figure 2, plot where points which end up in the Npnt points
#started from. NOTE -- this only includes pathways that start in the
#habitat! Note that the code is very similar to that for figure 1, but
#getting the inverse connectivity from Etranspose, and with nxTo
#becoming nxFrom and vice versa, and numTo become numFrom
figure(2);clf();style.use('ggplot')

#IMPORTANT, there is no gaurentee that the ordering of points in E and
#Etranspose are the same. In other words, there is no gaurentee that
#E.nxFrom[3]=Etranspose.nxTo[3]. So we must recalculate closepoints
#for Etransform. So now find the Npoint To points in the habitat
#closest to lonP,latP. It would be more efficient and accurate to use a
#plane-sailing formula, but this should be good enough, is simple, and
#requires no external libraries
lonTo,latTo=nxny2lonlat(Etranspose.nxTo,Etranspose.nyTo) #lon,lat of arrival (To) points
scaleDistSq=((lonTo-lonP)*cos(deg2rad(latTo)))**2+(latTo-latP)**2 #in weird units
closePoints=argsort(scaleDistSq)[:Npnts]


#plot the location of the habitat, using the lonTo,latTo points
#calculated above. One could also lonTo and latTo from the Etransform
plot(lonTo,latTo,'k,',markersize=5,zorder=10)
plot(lonTo[closePoints],latTo[closePoints],'r*',zorder=11)

#now, for each point larvae can arrive to, there are many points they
#could have come from. So for each lonTo and latTo, there are many
#lonFrom and latFrom points. So while Etranspose.nxTo[3] returns a
#single point, Etranspose.nxFrom[3] and Etranspose.nyFrom[3] will
#return a list of nx and ny coordinantes, and Etranspose.numFrom[3]
#will return an integer giving how many drifters came from each
#point. So a default dictionary will be used to count up how many
#particles go to each point. This double loop is certainly not the
#most efficient way to do this. But to make it a little faster, lets
#load the To variables into memory
EnxFrom=Etranspose.nxFrom[:]; EnyFrom=Etranspose.nyFrom[:]; EnumFrom=Etranspose.numFrom[:]
totalFrom=defaultdict(lambda :0)
for nPoint in closePoints:
    for pNum in zip(EnxFrom[nPoint],EnyFrom[nPoint],EnumFrom[nPoint]):
        p=pNum[:2]
        totalFrom[p]=totalFrom[p]+pNum[2]

nxFromVecTotal=[p[0] for p in totalFrom]
nyFromVecTotal=[p[1] for p in totalFrom]
numTotal=[totalFrom[p] for p in totalFrom]

#now make a scatter plot of where they end up
lonFromTotal,latFromTotal=nxny2lonlat(nxFromVecTotal,nyFromVecTotal)
scatter(lonFromTotal,latFromTotal,s=5,c=numTotal,vmin=0,zorder=12)
colorbar()

title('where particles which came to red stars started,\n'+
      'only including particles released from habitat')
draw()
show()
