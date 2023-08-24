from pylab import *
from numpy import *
import xarray as xr
import netCDF4 as nc
import time
import sklearn.neighbors as skn
import zarr
import copy

#specify where the files that will end up in EZfateFiles are located
dataDir='/Volumes/pringleExternal01/EZfateFiles/'

#load model grid from Mercator grid file
gridFile='ext-PSY4V3R1_mesh_zgr.nc'
gridData=nc.Dataset(dataDir+gridFile,'r')

#get grid data
nav_lon=gridData['nav_lon'][:]
nav_lat=gridData['nav_lat'][:]
depth=gridData['hdepw'][0,:]

nav_nx,nav_ny=meshgrid(arange(nav_lon.shape[1]),arange(nav_lon.shape[0])) #note, arrays indexed (ny,nx)


#make a function to return lat and lon of points that meet some depth
#criteria
def getPointsByDepth(minDepth,maxDepth):
    #ocean depths are positive

    depth=squeeze(gridData['hdepw'][0,:,:])
    depth[:,0]=nan; depth[:,-1]=nan; depth[0,:]=nan;depth[-1,:]=nan

    indx=logical_and(depth>=minDepth,depth<=maxDepth)

    return nav_lon[indx],nav_lat[indx]

#make a function to return grid indices of points that meet some depth
#criteria
def getGridPointsByDepth(minDepth,maxDepth):
    #ocean depths are positive

    depth=squeeze(gridData['hdepw'][0,:,:])
    depth[:,0]=nan; depth[:,-1]=nan; depth[0,:]=nan;depth[-1,:]=nan

    indx=logical_and(depth>=minDepth,depth<=maxDepth)

    return nav_nx[indx],nav_ny[indx]


#we want to find distance from land or ocean shallower than some depth
#in order to more easily define habitats. So we iterate over the depths
#of particle releases.
depthNameVec=['','_10m','_20m','_40m']
depthVec=[1,10,20,40]

#now work on each of these, one a at a time
if False: #so we can turn off this bit for debugging the code that writes the files
    for nD in range(len(depthVec)):
        thisDepth=depthVec[nD]
        depthName=depthNameVec[nD]


        #get land and sea points. Take into account that land in this
        #model often has a depth of 2.8m for some reason
        maxDland=max(2.8,thisDepth)
        print('Assuming the maximum depth of ground is',maxDland)

        #now get points for land
        minD=0; maxD=maxDland
        lonVecLand,latVecLand=getPointsByDepth(minD,maxD)
        nxVecLand,nyVecLand=getGridPointsByDepth(minD,maxD)

        #now zip all land points into an array of the lat and lon of land
        #points and make a BallTree from them with a haversine metric (which
        #requires them to converted to radians.
        #NOTE THAT THE HAVERSINE METRIC REQUIRES THE ORDER (LAT,LON) AND RADIANS
        landPoints=zeros((len(lonVecLand),2))
        landPoints[:,1]=lonVecLand
        landPoints[:,0]=latVecLand
        print('making land tree for depthName',depthName)
        landTree=skn.BallTree(radians(landPoints),metric='haversine')

        #these runs take a while, so enable each to be turned off so we can debug the other
        if True: #compute distance from each point to land in km
            #now make vectors of all model points, and find distance from each of
            #those to closest land point
            #NOTE THAT THE HAVERSINE METRIC REQUIRES THE ORDER (LAT,LON) AND RADIANS
            print('finding distance to land in km for depthName',depthName)
            tic=time.time()
            Re=6371.0
            allPoints=zeros((prod(nav_lon.shape),2))
            allPoints[:,1]=nav_lon.flatten()
            allPoints[:,0]=nav_lat.flatten()
            dist,i=landTree.query(radians(allPoints),k=1)
            dist=Re*dist.reshape(nav_lon.shape) #reshape and convert from radians to km
            print('   done in',time.time()-tic)

            #now save dist to a zarr file
            zarr.save(dataDir+'modelDistance%s.zip'%(depthName,),dist=dist)


        #now zip all land points into an array of the nx and ny of land points
        #and make a BallTree. Do (ny,nx) order to preserve sanity with respect
        #to work above.
        landPoints=zeros((len(nxVecLand),2))
        landPoints[:,1]=nxVecLand
        landPoints[:,0]=nyVecLand
        print('making land grid tree')
        landGridTree=skn.BallTree(landPoints)


        if True: #compute distance in grid-cells
            print('finding distance to land in grid points for depthName',depthName)
            tic=time.time()
            allPoints=zeros((prod(nav_nx.shape),2))
            allPoints[:,1]=nav_nx.flatten()
            allPoints[:,0]=nav_ny.flatten()
            dist,i=landGridTree.query(allPoints,k=1)
            gridDist=dist.reshape(nav_lon.shape) 
            print('   done in',time.time()-tic)

            #now save dist to a zarr file
            zarr.save(dataDir+'modelGridDistance%s.zip'%(depthName,),gridDist=gridDist)

        print(' ')

#now if true, make into a single netCDF file
if True:
    print('now write out data')
    #now save all grid information out as an xarray dataset which has lat
    #and lon as coordinates, and depth, dist* and gridDist* as variables. 

    #set up dataSet using my favorite resource
    #https://towardsdatascience.com/how-to-create-xarray-datasets-cf1859c95921
    attrs={'created_by':'add_distance_to_MercatorGridFile.py',
           'fromMercatorFile':gridFile}
    data_vars=dict(depth=(["ny","nx"],depth)
    #               dist=(["ny","nx"],dist),
    #               gridDist=(["ny","nx"],gridDist)
    )
    coords=dict(nav_lon=(["ny","nx"],nav_lon),
                nav_lat=(["ny","nx"],nav_lat)
    )

    #now add to data_vars the variables for each depth
    for depthName in depthNameVec:
        print('getting depth',depthName)
        gridDist=zarr.load(dataDir+'modelGridDistance%s.zip'%(depthName,))['gridDist']#.astype(float32)
        dist=zarr.load(dataDir+'modelDistance%s.zip'%(depthName,))['dist']#.astype(float32)

        data_vars['gridDist'+depthName]=(["ny","nx"],gridDist.copy())
        data_vars['dist'+depthName]=(["ny","nx"],dist.copy())

    print('writing data out')
    dataOut=xr.Dataset(attrs=attrs,data_vars=data_vars,coords=coords)
    comp = dict(zlib=True, complevel=5) #what compression
    encoding = {var: comp for var in dataOut.data_vars} #make dict all variables

    #for stupid reasons, writing dataOut straight to netcdf is slow,
    #because it is not set up for dask. So lets save it to zarr, load
    #that back, and then write it out. WTF?
    dataOut.to_zarr(dataDir+'jnk.zarr')
    dataOut=xr.open_zarr(dataDir+'jnk.zarr')
    
    #see https://github.com/pydata/xarray/issues/2912 for the reason for the .load()
    dataOut.to_netcdf(dataDir+'model_depth_and_distance.nc',encoding=encoding)
    print('done')
