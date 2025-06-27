from pylab import *
from numpy import *
import xarray as xr
from collections import defaultdict
import time
import dask
from collections import Counter
import zarr
import numcodecs
import os
import netCDF4 as nc
import sklearn.neighbors as skn
import pandas as pd
import shutil
import copy
from shapely.geometry import Point,MultiPoint
from shapely.geometry.polygon import Polygon
import getEZfateFromOSN


# this module has code to make multiple conectivity matrices from one
# or more trajectory files as filtered by depth and distance from land
# for the starting points, as well as the time the drifter started
# (either real time, or season).
#
# The locations are stored as integers which map to the model
# grid. The code to do is grid2int and int2grid, and they must be
# adjusted if the model changes.

#===============================================================================================
#define locations of file that define the model grid and Lagrangian connectivity
#from the EZfate project, as described in https://github.com/JamiePringle/EZfate
#the function getEZfateFromOSN.getFileFromOSN() takes as an argument a pathway to the
#EZfate data on the open storage network S3 bucket, and returns the path to a local
#file that contains the same data. Details as to where and how this done, and where
#the data is stored, can be found in the getEZfateFromOSN module.

#the location of the mask grid for the model run. This defines the
#locations of all the model data files
maskFile='EZfateData/EZfateFiles/ext-PSY4V3R1_mesh_zgr.nc'
maskFile=getEZfateFromOSN.getFileFromOSN(maskFile)

#where are dispersal matrices
matDir='EZfateData/communityConnectivityMatrices/'
matDir=getEZfateFromOSN.getFileFromOSN(matDir)

#===============================================================================================
#define mapping from grid indices nx,ny to a single number
#the following take arrays as arguements
pointType=uint32
def grid2int(nx,ny):
    '''
    grid2int(nx,ny):

    nx and ny are arrays representing model grid points

    convert nx,ny to a single integer
    '''
    out=empty(nx.shape,dtype=pointType)
    out[:]=pointType(nx)*10000+pointType(ny)
    return out

def int2grid(points):
    '''int2grid(points): convert from a single integer array to mdoel
    grid points to nx and ny; assumes points is an integer or array

    '''
    nx=points//10000
    ny=points-nx*10000
    return nx,ny

#===============================================================================================
# routines to get
#    -depth of water of each point getDepthDict()
#    -distance from shore of each point in kilometers getDistanceDict()
#    -distance from shore of each point in grid spacing getGridDistanceDict()

def getDepthDict(nxVec,nyVec):
    '''
    getDepthDict(nxVec,nyVec):
    #this code takes a vector grid indices nxVec and nyVec and returns
    #a dictionary that maps from (nx,ny) to the depth at that grid
    #point
    '''
    gridData=nc.Dataset(maskFile)
    depth=gridData['hdepw'][0,:,:]
    depthDict={(p[0],p[1]):depth[p[1],p[0]] for p in zip(nxVec,nyVec) }
    return depthDict

def getDistanceDict(nxVec,nyVec,landThresh=2.1):
    '''getDistanceDict(nxVec,nyVec): this code takes a vector grid indices
    nxVec and nyVec and returns a dictionary that maps from (nx,ny) to
    distance from land in km. Land is defined as any distance less
    than landThresh, which defaults to 2.1 as appropriate for the
    Mercator 1/12th global run

    '''
    gridData=nc.Dataset(maskFile)
    depth=gridData['hdepw'][0,:,:]
    nav_lon=gridData['nav_lon'][:]
    nav_lat=gridData['nav_lat'][:]

    #now make coast grid
    #landThresh=2.1
    print('Assuming land is depth<',landThresh)
    landGrid=depth<landThresh
    coastGrid=logical_and(False,landGrid) #start coast as all false
    coastGrid[1:-1,1:-1]=logical_and(landGrid[1:-1,1:-1], #coast is a point that is land next to a point that is water
                                     logical_not(functools.reduce(logical_and,(landGrid[:-2,2:],
                                                                               landGrid[1:-1,2:],
                                                                               landGrid[2:,2:],
                                                                               landGrid[:-2,1:-1],
                                                                               landGrid[2:,1:-1],
                                                                               landGrid[:-2,:-2],
                                                                               landGrid[1:-1,:-2],
                                                                               landGrid[2:,:-2]
                                     ))))

    coastLon=nav_lon[coastGrid]
    coastLat=nav_lat[coastGrid]
    coastPoints=zeros((len(coastLon),2))
    coastPoints[:,1]=coastLon
    coastPoints[:,0]=coastLat
    print('making coast point tree')
    coastTree=skn.BallTree(radians(coastPoints),metric='haversine')
    allPointsLatLon=array([(nav_lat[p[0],p[1]],nav_lon[p[0],p[1]]) for p in zip(nyVec,nxVec)])
    dist,jnk=coastTree.query(radians(allPointsLatLon),k=1,return_distance=True)
    Re=6371.0 #in km
    dist=Re*dist
    distanceDict={(p[0],p[1]):p[2][0] for p in zip(nxVec,nyVec,dist) }
    
    return distanceDict

def getGridDistanceDict(nxVec,nyVec,landThresh=2.1):
    '''getGridDistanceDict(nxVec,nyVec): this code takes a vector of grid
    indices nxVec and nyVec and returns a dictionary that maps from
    (nx,ny) to distance from land in grid spacings. Land is defined as
    any distance less than landThresh, which defaults to 2.1 as
    appropriate for the Mercator 1/12th global run

    '''
    gridData=nc.Dataset(maskFile)
    depth=gridData['hdepw'][0,:,:]

    nav_nx,nav_ny=meshgrid(arange(depth.shape[1]),arange(depth.shape[0]))

    #now make coast grid
    #landThresh=2.1
    print('Assuming land is depth<',landThresh)
    landGrid=depth<landThresh
    coastGrid=logical_and(False,landGrid) #start coast as all false
    coastGrid[1:-1,1:-1]=logical_and(landGrid[1:-1,1:-1], #coast is a point that is land next to a point that is water
                                     logical_not(functools.reduce(logical_and,(landGrid[:-2,2:],
                                                                               landGrid[1:-1,2:],
                                                                               landGrid[2:,2:],
                                                                               landGrid[:-2,1:-1],
                                                                               landGrid[2:,1:-1],
                                                                               landGrid[:-2,:-2],
                                                                               landGrid[1:-1,:-2],
                                                                               landGrid[2:,:-2]
                                     ))))

    coastNx=nav_nx[coastGrid]
    coastNy=nav_ny[coastGrid]
    coastPoints=zeros((len(coastNx),2))
    coastPoints[:,1]=coastNx
    coastPoints[:,0]=coastNy
    print('making coast point tree')
    coastTree=skn.BallTree(coastPoints)
    print('   done with tree')
    allPointsNyNx=array([(nav_ny[p[0],p[1]],nav_nx[p[0],p[1]]) for p in zip(nyVec,nxVec)])
    dist,jnk=coastTree.query(allPointsNyNx,k=1,return_distance=True)
    gridDistanceDict={(p[0],p[1]):p[2][0] for p in zip(nxVec,nyVec,dist) }
    
    return gridDistanceDict



#===============================================================================================
# make the simplest connectivity matrix. This version of the code save
# the connectivity matrix as a zarr file. Each row represents an
# origin point, and holds all the destination points. The coordinates
# are model grid points (nx,ny) in the drifter data, and are stored as
# type coordType. coordtype is now uint16.  The variables are
#
# nxFrom: grid index of origin of drifter
# nyFrom: grid index of origin of drifter
# nxTo: grid indices of ending points of drifters, stored as ragged array
# nyTo: grid indices of ending points of drifters, stored as ragged array
# numTo: number of points that went to each ending point, stored as ragged array

#type of model coordinate.
coordType=uint16

# the type of the number of drifters that make a particular
# trip. Warning, danger of overflow...
typeOfSum=uint16
def typeOfColumn():
    return defaultdict(typeOfSum)

#===============================================================================================
def makeConnectFromFileWithinHabitat(data,ageIndex,goodHabitat,dirOut,minDay=-1, maxDay=500,
                                     onlySettleInside=True):
    '''makeConnectFromFile(data,ageIndex,dirOut,minDay=-1, maxDay=500, onlySettleInside=True):

    This function makes the connectivity matrix for all points
    starting in the supplied data that start and end within the points
    defined by goodHabitat. If multiple ages are provided in ageIndex,
    the first one in the habitat defines the connectivity.

    It assumes that the input drifter data is set up so each coloumn
    has the same time since release, and will limit the data used to
    drifters starting on year days between minDay and maxDay. If
    ageIndx includes multiple times, it will include the drifter
    multiple times, once for each age in ageIndex.

    Every point included in goodHabitat will exist in the output
    connectivity matrix -- however, the nxTo,nyTo and numTo values
    may be empty.

    if onlySettleInside=True (default), then particles will only
    settle in points in goodHabitat.  If False, They can settle
    anywhere. If False, then ageIndex should only have one true value,
    because particles will always in first possible time. 

    data: the xarray/zarr dataset that defines the drifter data
    (as made by 03_addGridPoints.py)

    ageIndex: a boolean index of the ages to include in the
    connectivity matrix whose size is the same data.dims['obs'] and
    whose values are true for the ages to be included. Connectivity is
    happens at the smallest time in ageIndex when the drifter is in
    the habitat.

    goodHabitat is a set whose keys are (nx,ny). If a point is in the
    set, they are included in the connectivity matrix.

    dirOut: the directory to which the connectivity matrix will be
    written as a zarr array.

    minDay=-1, maxDay=500: only drifters which start between year-day
    minDay and maxDay are included in the connectivity matrix; the
    default values will include the entire year.

    does not return anything.

    '''
    #now lets pretend we can get the data into memory all at once, and
    #load it from disk
    tic=time.time()

    #to get in one chunk, get all destination ages and the start
    #location (obs=0) all at once
    ageIndexAll=ageIndex.copy()
    ageIndexAll[0]=True

    print('getting all date data')
    if False:
        #only works with xarray which keeps track of cf date conventions 
        startDay=pd.DatetimeIndex(data['time'][:,0]).dayofyear #yearday of start of drifter
    else:
        startDay=data['time'][:,0]/8.64e4+1.0
        print('MAKING DANGEROUS ASSUMPTION THAT FLOAT RUN STARTED ON FIRST DAY OF THE YEAR')
        
    indxStartToKeep=logical_and(startDay>=minDay,startDay<=maxDay)

    #find number of starting times, to write into output file
    numberOfStarts=len(unique(startDay[indxStartToKeep]))

    print('getting all position data')
    pointsAll=grid2int(data['nx'].oindex[indxStartToKeep,ageIndexAll],
                         data['ny'].oindex[indxStartToKeep,ageIndexAll])
    pointsStart=pointsAll[:,0]
    pointsLater=pointsAll[:,1:]
    del pointsAll
    print('   done in',time.time()-tic,'with',len(pointsStart),'drifters out of',len(startDay),'in file')

    #Now, some drifters will have gone through a boundary and
    #died. Their nx and ny will be maxInt, or about 65535 and after a
    #grid2int(), the value will be 655415535. I am nervous this will
    #introduce some archetecture dependence... but check for too big
    #values of pointsLater, and delete them. I wish int's had NaNs
    indxGoodPoints=(amax(pointsLater,axis=1)<650000000)
    pointsLater=pointsLater[indxGoodPoints,:]
    pointsStart=pointsStart[indxGoodPoints]

    #now I want to find all drifters that enter good habitat in one of
    #the times in pointsLater, and also which point they first
    #entered. Lets map goodHabitat into a single integer
    goodPoints=set()
    goodPoints.update([grid2int(array([p[0]]),array([p[1]]))[0] for p in goodHabitat])

    #first trim all points that did not start in habitat.
    indx=array([(p in goodPoints) for p in pointsStart])
    pointsStart=pointsStart[indx]
    pointsLater=pointsLater[indx,:]
    print('done finding those that started in habitat, n=',len(pointsStart))

    #now march from last to first time, and see which of pointsLater
    #were first to be in goodPoints (or anywhere, if onlySettleInside=False
    inHabitat=pointsStart<0 #should be all false
    firstGood=0*pointsStart+9999 #initialize this to be the right size.
    for n in flip(arange(sum(ageIndex))):
        print('   looking at time i=',n)
        if onlySettleInside:
            isIn=array([(p in goodPoints) for p in pointsLater[:,n]])
        else:
            isIn=array([True for p in pointsLater[:,n]]) #all True
        inHabitat[isIn]=True
        firstGood[isIn]=n

    print('done finding those which end in habitat, n=',sum(isIn))

    #breakpoint()
    
    #now provide all points where particles start and end in habitat,
    #giving the first time they enter the habitat.
    pointsStartDontEnd=pointsStart[logical_not(isIn)] # this will be used to make empty bits of connectivity
    pointsStart=pointsStart[isIn]
    pointsEnd=pointsLater[isIn,:] #keep only those in
    firstGood=firstGood[isIn] #keep only those that are in somewhere
    #pointsEnd=pointsEnd[:,firstGood] #keep the first one to be in
    pointsEnd=array([pointsEnd[n,firstGood[n]] for n in range(len(firstGood))])

    print('Shape of pointsStart and pointsEnd are',pointsStart.shape,pointsEnd.shape)

    #now iterate over ages to get
    tic=time.time()
    #breakpoint()
    #try counters, with key being (from,to)
    #about 25s/age but now in wrong shape.

    #make connectivity matrix.
    E=Counter()
    ticInner=time.time()
    E.update(zip(pointsStart,pointsEnd))
    print('done making E')

    #now we want to make a structure that Eout[nFrom][nTo] gives
    #the number of points that went from nFrom to nTo. This is
    #done in a tight little list comprehension which I need to
    #debug carefully. Warning, read about Counter.update carefully.  It is
    #funny, and requires a mapping to work the way you might
    #expect
    ticInner=time.time()
    Eout=defaultdict(Counter)
    [Eout[x[0]].update({x[1]:E[x]}) for x in E ] #if x[0]==31232452]
    print('   done making Eout in',time.time()-ticInner)

    #breakpoint()
        
    print('   done organizing data in ',time.time()-tic,'now make zarr array')
    tic=time.time()
    
    #start making zarr array in dirOut. 

    
    nOut=len(Eout) #This is the number of points that connected to something
    nTotal=len(goodHabitat) #This is the number of points in the starting habitat

    #find all starting points from goodHabitat that are not in output
    #because non of their offspring settled in habitat
    notInPoints=goodPoints-Eout.keys()
    assert nTotal-nOut==len(notInPoints),'Sanity check failed, think'

    #create zarr dataset
    chunkSize=int(1e6) #should experiment with this
    #store=zarr.DirectoryStore(dirOut) #write as directory
    store=zarr.ZipStore(dirOut) #write as zip 
    root=zarr.group(store=store)
    
    nxFrom=root.empty(shape=(nTotal,),name='nxFrom',dtype=coordType,chunks=chunkSize)
    nyFrom=root.empty(shape=(nTotal,),name='nyFrom',dtype=coordType,chunks=chunkSize)

    nxTo=root.empty(shape=(nTotal,),name='nxTo',dtype=object,
                    chunks=chunkSize,object_codec=numcodecs.VLenArray(coordType))
    nyTo=root.empty(shape=(nTotal,),name='nyTo',dtype=object,
                    chunks=chunkSize,object_codec=numcodecs.VLenArray(coordType))
    numTo=root.empty(shape=(nTotal,),name='numTo',dtype=object,
                     chunks=chunkSize,object_codec=numcodecs.VLenArray(typeOfSum))

    #write number of particles which are released as an attribute
    root.attrs['numberOfStartingTimes']=int(numberOfStarts)
    
    #The following code is slow, and assumes that everytime you read a
    #dictionary it gives you the answers in the same order, which is
    #only true for recent versions of python.

    #do nxFrom and nyFrom
    print('writing out nxFrom,nyFrom')
    nx,ny=int2grid(array(list(Eout)))
    nxFrom[:nOut]=nx
    nyFrom[:nOut]=ny

    nxEmpty,nyEmpty=int2grid(array(list(notInPoints)))
    nxFrom[nOut:]=nxEmpty
    nyFrom[nOut:]=nyEmpty
    print('   done')

    #breakpoint()
    
    #First, write out the points which started and connected to something
    nIn=0 #counter for number of records written out
    nxToNow=[] #to store data before writing
    nyToNow=[]  #to store data before writing
    numToNow=[]  #to store data before writing
    for key in Eout:
        nIn+=1

        nItem=len(Eout[key])
        nxToNowTemp,nyToNowTemp=int2grid(array(list(Eout[key])))
        nxToNow.append(nxToNowTemp)
        nyToNow.append(nyToNowTemp)
        numToNow.append(array([Eout[key][e] for e in Eout[key]]))

        if amax(nyToNowTemp)>5000 : #THIS IS FOR DEBUGGING ONLY, SHOULD NOT BE IN PRODUCTION CODE
            breakpoint()
            

        if (len(nxToNow)==chunkSize) or (nIn==nOut):
            #write out list
            #cOut=min(nIn+chunkSize,nOut)
            #breakpoint()
            cOut=len(nxToNow)
            if nIn==nOut:
                #sigh, wrestling with how array works. When there
                #is just one array in list, it is a struggle to make
                #it in an array of type object with just one element
                jnk=zeros((cOut,),dtype=object); jnk[:]=nxToNow; nxTo[nIn-cOut:nIn]=jnk
                jnk=zeros((cOut,),dtype=object); jnk[:]=nyToNow; nyTo[nIn-cOut:nIn]=jnk
                jnk=zeros((cOut,),dtype=object); jnk[:]=numToNow; numTo[nIn-cOut:nIn]=jnk
            else:
                #can't insert list into an array
                nxTo[nIn-cOut:nIn]=array(nxToNow,dtype=object)
                nyTo[nIn-cOut:nIn]=array(nyToNow,dtype=object)
                numTo[nIn-cOut:nIn]=array(numToNow,dtype=object)
            #if remainder(nIn,100)==0:
            #    print(nIn/nOut*100,'percent done in nxTo,nyto and numTo')
            nxToNow=[] #to store data before writing
            nyToNow=[]  #to store data before writing
            numToNow=[]  #to store data before writing


    #close store. Important for ZipStore
    store.close()

    print('Done with writing in ',time.time()-tic)

    return None


#===============================================================================================
# adapt the connectivity routine above to make a connectivity matrix
# that only includes points starting within a certain depth and
# distance from shore range, and in a certain time of year range

#===============================================================================================
# cull from connectivity matrices points that start or stop at
# locations that do not meet some criterion.

def isGood(theDict,theMin,theMax,thePoint):
    '''returns True if thePoint is in theDict AND theMin<=theDict[thePoint]<=theMax'''
    #note that this works because of how python short-circuits evaluation of 'if'
    return (thePoint in theDict) and (theMin<=theDict[thePoint]) and (theDict[thePoint]<=theMax)

#@profile
def trimConnectivity(connectivityIn,theDict,theMin,theMax,loadConnectToMem=True,keepAllTo=False):
    '''trimConnectivity(connectivityIn,theDict,theMin,theMax) takes as input

    connectivityIn, an open zarr file that specifies the connectivity
    matrix, as created by makeConnectFromFile()

    theDict, a dictionary with key (nx,ny) that returns a number as a
    value (e.g. as from depthDict(), where the dictionary returns the
    depth at the location (nx,ny).

    theMin and the theMax: each starting and ending point in
    connectivityIn (defined by their (nx,ny)) is only included in the
    new connectivity matrix if (nx,ny) is in theDict AND
    theMin<=theDict(nx,ny)<=theMax

    loadConnectToMem=True: If this is true, load connectivity into
    memory -- much faster, but much more memory. If false, work on disk, slower but no memory ussage.

    keepAllTo=False: if False, also trim *To points that are meet the triming criteria. if True, don't
    remove *To points

    returns connectivityOut, an in memory zarr file that is the
    reduced connectivity matrix that satisfies the criterion above.

    '''

    #first, make in-memory connectivityOut zarr data store and fill
    #with appropriate (nxFrom,nyFrom)
    nxFromIn=connectivityIn['nxFrom'][:]
    nyFromIn=connectivityIn['nyFrom'][:]

    toKeep=[isGood(theDict,theMin,theMax,p) for p in zip(nxFromIn,nyFromIn)]

    #make output as an in-memory zarr
    nOut=sum(toKeep)

    #create zarr dataset
    chunkSize=int(1e6) #should experiment with this
    store=zarr.MemoryStore()
    root=zarr.group(store=store)
    
    nxFrom=root.empty(shape=(nOut,),name='nxFrom',dtype=coordType,chunks=chunkSize)
    nyFrom=root.empty(shape=(nOut,),name='nyFrom',dtype=coordType,chunks=chunkSize)
    if ('numlaunched' in connectivityIn) or ('numLaunched' in connectivityIn): #zarr lowercases everything
        print('   adding numLaunched variable')
        numLaunched=root.empty(shape=(nOut,),name='numLaunched',dtype='i',chunks=chunkSize)

    nxTo=root.empty(shape=(nOut,),name='nxTo',dtype=object,chunks=chunkSize,
                    object_codec=numcodecs.VLenArray(coordType))
    nyTo=root.empty(shape=(nOut,),name='nyTo',dtype=object,chunks=chunkSize,
                    object_codec=numcodecs.VLenArray(coordType))
    numTo=root.empty(shape=(nOut,),name='numTo',dtype=object,chunks=chunkSize,
                     object_codec=numcodecs.VLenArray(typeOfSum))

    #write out nxFrom,nyFrom
    nxFrom[:]=nxFromIn[toKeep]
    nyFrom[:]=nyFromIn[toKeep]
    if ('numlaunched' in connectivityIn) or ('numLaunched' in connectivityIn): #zarr lowercases things
        print('  adding numLaunched data')
        numLaunchedIn=connectivityIn['numLaunched'][:] #again, zarr lowercases things? or not? or macOS
        numLaunched[:]=numLaunchedIn[toKeep] 

    #now we need to update nxTo,nyTo and numTo
    print('starting to fill lists')
    tic=time.time()
    sourceIndx=arange(len(nxFromIn))[toKeep]
    nxToNew=empty((nOut,),dtype=object)
    nyToNew=empty((nOut,),dtype=object)
    numToNew=empty((nOut,),dtype=object)

    #if we can fit the connectivity array in memory, it should be much
    #faster... if we can't, just leave as link to zarr array
    if not loadConnectToMem:
        #slow, low memory usage
        print('In trimConnectivity, using slow, low memory ussage path')
        nxToIn=connectivityIn['nxTo']
        nyToIn=connectivityIn['nyTo']
        numToIn=connectivityIn['numTo']
    else:
        #fast, high memory usage
        print('In trimConnectivity, Loading connectivity matrix into memory; if fails, make loadConnectToMem=False')
        nxToIn=connectivityIn['nxTo'][:]
        nyToIn=connectivityIn['nyTo'][:]
        numToIn=connectivityIn['numTo'][:]
        print('   done loading into memory')

    #breakpoint()
    
    for n in range(nOut):
        nIn=sourceIndx[n]
        nxToOld=nxToIn[nIn]
        nyToOld=nyToIn[nIn]
        numToOld=numToIn[nIn]
        if keepAllTo:
            nxToNew[n]=nxToOld
            nyToNew[n]=nyToOld
            numToNew[n]=numToOld
        else:
            toKeep=[isGood(theDict,theMin,theMax,p) for p in zip(nxToOld,nyToOld)]
            nxToNew[n]=nxToOld[toKeep]
            nyToNew[n]=nyToOld[toKeep]
            numToNew[n]=numToOld[toKeep]
    print('done filling lists in ',time.time()-tic)
    nxTo[:]=nxToNew 
    nyTo[:]=nyToNew
    numTo[:]=numToNew
    print('done with To',time.time()-tic,'cummalative')

    #pass the numberOfStartingTimes info through to the output
    root.attrs['numberOfStartingTimes']=connectivityIn.attrs['numberOfStartingTimes']
    
    return root

#===============================================================================================
# cull from connectivity matrices points that start or stop at
# locations that are inside or outside of some polygon

# def isGood_forPoly(theDict,thePoint):
#     '''returns True if thePoint is in theDict AND theDict[thePoint] is true'''
#     #note that this works because of how python short-circuits evaluation of 'if'
#     return (thePoint in theDict) and theDict[thePoint]

#@profile
#print('TURNIPS!, TURNIPS')
def trimConnectivity_byPoly(connectivityIn,thePolygon,keepInPoly=True,loadConnectToMem=True,keepAllTo=False):
    '''trimConnectivity(connectivityIn,theDict,theMin,theMax) takes as input

    connectivityIn, an open zarr file that specifies the connectivity
    matrix, as created by makeConnectFromFile()

    thePolygon, a list of (lon,lat) points that define a polygon.

    keepInPoly=True, if this is true, keep points in polygon, if false, keep points outside polygon

    loadConnectToMem=True: If this is true, load connectivity into
    memory -- much faster, but much more memory. If false, work on disk, slower but no memory ussage.

    keepAllTo=False: if False, also trim *To points that are meet the triming criteria. if True, don't
    remove *To points

    returns connectivityOut, an in memory zarr file that is the
    reduced connectivity matrix that satisfies the criterion above.

    '''

    #first, make in-memory connectivityOut zarr data store and fill
    #with appropriate (nxFrom,nyFrom)
    nxFromIn=connectivityIn['nxFrom'][:]
    nyFromIn=connectivityIn['nyFrom'][:]

    #get latitude and longitude matrices, and make vectors of latFromIn and lonFromIn
    with nc.Dataset(maskFile,'r') as gridData:
        lonMat=gridData['nav_lon'][:]
        latMat=gridData['nav_lat'][:]
        lonHabVec=array([lonMat[p] for p in zip(nyFromIn,nxFromIn)])
        latHabVec=array([latMat[p] for p in zip(nyFromIn,nxFromIn)])

    #now make dictionary with key (nx,ny) and a boolean if in polygon
    allPoints=MultiPoint([a for a in zip(lonHabVec,latHabVec)])
    regionPolyShapely=Polygon(thePolygon)
    toKeep=[p.within(regionPolyShapely) for p in allPoints.geoms] #there has to be a better way

    #plot(*regionPolyShapely.exterior.xy)
    #
    # jnkx=[p.x for p in allPoints.geoms]
    # jnky=[p.y for p in allPoints.geoms]
    # plot(jnkx,jnky,'k*')
        
                        
    #if keepInPoly is false, then we want to keep the points outside
    #of the polygon, so must reverse toKeep
    if not keepInPoly:
        toKeep=logical_not(toKeep)

    #make dictionary of points to keep with key (nx,ny)
    toKeepDict={}
    jnk=[(p,True) for p in zip(nxFromIn[toKeep],nyFromIn[toKeep])]
    toKeepDict.update(jnk)

    #make output as an in-memory zarr
    nOut=sum(toKeep)

    #create zarr dataset
    chunkSize=int(1e6) #should experiment with this
    store=zarr.MemoryStore()
    root=zarr.group(store=store)

    #print('   LOST IN TURNIPS'); breakpoint()
    
    nxFrom=root.empty(shape=(nOut,),name='nxFrom',dtype=coordType,chunks=chunkSize)
    nyFrom=root.empty(shape=(nOut,),name='nyFrom',dtype=coordType,chunks=chunkSize)
    if ('numlaunched' in connectivityIn) or ('numLaunched' in connectivityIn): #zarr lowercases everything
        print('   adding numLaunched variable')
        numLaunched=root.empty(shape=(nOut,),name='numLaunched',dtype='i',chunks=chunkSize)

    nxTo=root.empty(shape=(nOut,),name='nxTo',dtype=object,chunks=chunkSize,
                    object_codec=numcodecs.VLenArray(coordType))
    nyTo=root.empty(shape=(nOut,),name='nyTo',dtype=object,chunks=chunkSize,
                    object_codec=numcodecs.VLenArray(coordType))
    numTo=root.empty(shape=(nOut,),name='numTo',dtype=object,chunks=chunkSize,
                     object_codec=numcodecs.VLenArray(typeOfSum))

    #write out nxFrom,nyFrom
    nxFrom[:]=nxFromIn[toKeep]
    nyFrom[:]=nyFromIn[toKeep]
    if ('numlaunched' in connectivityIn) or ('numLaunched' in connectivityIn): #zarr lowercases things
        print('  adding numLaunched data')
        #is the case insensitivity of macOS screwing me up? Will this work on my laptop?
        numLaunchedIn=connectivityIn['numLaunched'][:] #again, zarr lowercases things
        numLaunched[:]=numLaunchedIn[toKeep] 


    #now we need to update nxTo,nyTo and numTo
    print('starting to fill lists')
    tic=time.time()
    sourceIndx=arange(len(nxFromIn))[toKeep]
    nxToNew=empty((nOut,),dtype=object)
    nyToNew=empty((nOut,),dtype=object)
    numToNew=empty((nOut,),dtype=object)

    #if we can fit the connectivity array in memory, it should be much
    #faster... if we can't, just leave as link to zarr array
    if not loadConnectToMem:
        #slow, low memory usage
        print('In trimConnectivity, using slow, low memory ussage path')
        nxToIn=connectivityIn['nxTo']
        nyToIn=connectivityIn['nyTo']
        numToIn=connectivityIn['numTo']
    else:
        #fast, high memory usage
        print('In trimConnectivity, Loading connectivity matrix into memory; if fails, make loadConnectToMem=False')
        nxToIn=connectivityIn['nxTo'][:]
        nyToIn=connectivityIn['nyTo'][:]
        numToIn=connectivityIn['numTo'][:]
        print('   done loading into memory')

    #breakpoint()
    
    for n in range(nOut):
        nIn=sourceIndx[n]
        nxToOld=nxToIn[nIn]
        nyToOld=nyToIn[nIn]
        numToOld=numToIn[nIn]
        if keepAllTo:
            nxToNew[n]=nxToOld
            nyToNew[n]=nyToOld
            numToNew[n]=numToOld
        else:
            toKeep=[(p in toKeepDict)  for p in zip(nxToOld,nyToOld)]
            nxToNew[n]=nxToOld[toKeep]
            nyToNew[n]=nyToOld[toKeep]
            numToNew[n]=numToOld[toKeep]
    print('done filling lists in ',time.time()-tic)
    nxTo[:]=nxToNew 
    nyTo[:]=nyToNew
    numTo[:]=numToNew
    print('done with To',time.time()-tic,'cummalative')

    #pass the numberOfStartingTimes info through to the output
    root.attrs['numberOfStartingTimes']=connectivityIn.attrs['numberOfStartingTimes']
    
    return root


#===============================================================================================
# Make code to combine connectivity matrices

def makeEmptyConnectivity(store):
    '''makeEmptyConnectivity(store)

    Take a zarr store (as created by zarr.MemoryStore() or
    zarr.DirectoryStore() or the like) and create an empty
    connectivity matrix in it. This can then be used by
    combineConnectivity to combine multiple matrices together.

    returns the root object to the connectivity matrix.

    '''
    #create zarr dataset
    chunkSize=int(1e6) #should experiment with this
    root=zarr.group(store=store)
    
    nOut=0 #empty!
    nxFrom=root.empty(shape=(nOut,),name='nxFrom',dtype=coordType,chunks=chunkSize)
    nyFrom=root.empty(shape=(nOut,),name='nyFrom',dtype=coordType,chunks=chunkSize)

    nxTo=root.empty(shape=(nOut,),name='nxTo',dtype=object,chunks=chunkSize,
                    object_codec=numcodecs.VLenArray(coordType))
    nyTo=root.empty(shape=(nOut,),name='nyTo',dtype=object,chunks=chunkSize,
                    object_codec=numcodecs.VLenArray(coordType))
    numTo=root.empty(shape=(nOut,),name='numTo',dtype=object,chunks=chunkSize,
                     object_codec=numcodecs.VLenArray(typeOfSum))

    return root

#@profile
def combineConnectivity(A,B,runFast=True):
    '''combineConnectivity(A,B,runFast=True)

    connectivity matrix B is added to connectivity in A. To combine
    many matrices, start with A as an empty matrix defined by
    makeEmptyConnectivity(), and keep adding to it.

    A and B are the root object to a zarr structure

    if runFast=True, then load some of B into memory, so that it is
    faster (but uses more memory)

    Does not return anything.

    '''

    #find total number of starting times in A and B, and update
    #numberOfStartingTimes in A. Need to deal with fact that A might not
    #have that attribute yet, if it is an empty dataset
    if 'numberOfStartingTimes' in A.attrs:
        numA=A.attrs['numberOfStartingTimes']
    else:
        numA=0
        
    numB=B.attrs['numberOfStartingTimes']
    A.attrs['numberOfStartingTimes']=numA+numB
        
    #make sets of starting points in A and B
    inA=set([(p[0],p[1]) for p in zip(A['nxFrom'][:],A['nyFrom'][:])])
    inB=set([(p[0],p[1]) for p in zip(B['nxFrom'][:],B['nyFrom'][:])])

    #make a dictionary that maps from (nxFrom,nyFrom) in B to row in B
    nxFromB=B['nxFrom'][:]
    nyFromB=B['nyFrom'][:]
    whereInB=dict([((nxFromB[n],nyFromB[n]),n) for n in range(len(nxFromB))])

    #find all points in A and B, and all points in B but not A
    inAandB=inA.intersection(inB)
    inBnotA=inB.difference(inA)

    if runFast:
        nxToAllB=B['nxTo'][:]
        nyToAllB=B['nyTo'][:]
        numToAllB=B['numTo'][:]

    #The basic logic is simple, but the implementation here is perhaps
    #slow.
    #
    #First loop over all points in A, by chunk, and if the starting
    #points are in A and B ((nxFrom,nyFrom) in inAandB) then add to
    #existing nxTo,nyTo and numTo. Note their are two cases to worry
    #about -- (nxTo,nyTo) in both A and B, and (nxTo,nyTo) just in
    #B. The case where a point is in A but not B should be fine and
    #automatically delt with in the follow loop for the in A and B
    #case.
    #
    #All this code needs to do is update nxTo,nyTo and numTo
    tic=time.time()
    print('Starting to combine inAandB')
    lenA=A['nxFrom'].shape[0]
    chunkSize=A['nxFrom'].chunks[0]
    startThisChunk=0
    while startThisChunk <lenA:
        thisChunkSize=min(chunkSize,lenA-startThisChunk) #on last loop, should == lenA 

        #get a chunk. Do a deep copy for the ones that will be
        #altered. Why? For reasons I don't understand, the object
        #arrays in *To arrays are unwriteable with the deep copy. It
        #is ok to do this, because below I copy the update chunks back
        #to A
        nxFrom=A['nxFrom'][startThisChunk:startThisChunk+thisChunkSize]
        nyFrom=A['nyFrom'][startThisChunk:startThisChunk+thisChunkSize]

        nxTo=copy.deepcopy(A['nxTo'][startThisChunk:startThisChunk+thisChunkSize])
        nyTo=copy.deepcopy(A['nyTo'][startThisChunk:startThisChunk+thisChunkSize])
        numTo=copy.deepcopy(A['numTo'][startThisChunk:startThisChunk+thisChunkSize])

        #update info in this chunk
        for nA in range(len(nxFrom)): #loop over chunk, nA is a row in the chunk
            if (nxFrom[nA],nyFrom[nA]) in inAandB: #is (nxFrom,nyFrom) in A and B

                #to check if To point from B in A already for this
                #(nxFrom,nyFrom), and to find where that
                #(nxFrom,nyFrom) are, create dictionary whose key is a
                #single (nxFrom[],nyFrom[]) in A and whose value is the
                #location of that point in the nxFrom/nyFrom/numFrom
                #vector for a row in A
                inThisTo=dict([((p[0],p[1]),p[2]) for p in zip(nxTo[nA],nyTo[nA],range(len(nxTo[nA])))]) 

                #breakpoint()
                
                #to store points that need to be appended to nxTo, nyTo and numTo of A
                nxToNew=[]; nyToNew=[]; numToNew=[]

                #loop over *To in B. Note that here it is gaurenteed
                #that *From points are in both A and B
                bRow=whereInB[(nxFrom[nA],nyFrom[nA])]
                if runFast:
                    nxToB=nxToAllB[bRow]
                    nyToB=nyToAllB[bRow]
                    numToB=numToAllB[bRow]                    
                else:
                    nxToB=B['nxTo'][bRow]
                    nyToB=B['nyTo'][bRow]
                    numToB=B['numTo'][bRow]

                for nB in range(len(nxToB)): #nB is, confusingly, a location in nxToB/nyToB/numToB
                    if (nxToB[nB],nyToB[nB]) in inThisTo:
                        #then the (nxTo,nyTo) point is in both A and B
                        #for this (nxFrom,nyFrom), and so we just need
                        #to update the numTo
                        whichOne=inThisTo[(nxToB[nB],nyToB[nB])]
                        numTo[nA][whichOne]+=numToB[nB]
                    else:
                        #then the (nxTo,nyTo) point is only in B for
                        #this (nxFrom,nyFrom), and so we need to add
                        #new items to nxTo/nyTo/numTo in A
                        nxToNew.append(nxToB[nB])
                        nyToNew.append(nyToB[nB])
                        numToNew.append(numToB[nB])

                #now, we need to append nxToNew,nyToNew and numToNew to nxTo[n]/nyTo[n]/numTo[n]
                nxTo[nA]=concatenate((nxTo[nA],array(nxToNew)))#,casting='no')
                nyTo[nA]=concatenate((nyTo[nA],array(nyToNew)))#,casting='no')
                numTo[nA]=concatenate((numTo[nA],array(numToNew)))#,casting='no')
                #breakpoint()

        #now write the chunk back to A, but only the *To's, since they
        #are the only ones to have changed
        A['nxTo'][startThisChunk:startThisChunk+thisChunkSize]=nxTo
        A['nyTo'][startThisChunk:startThisChunk+thisChunkSize]=nyTo
        A['numTo'][startThisChunk:startThisChunk+thisChunkSize]=numTo

        #now update startThisChunk
        #print('Done with',startThisChunk,'of',lenA)
        startThisChunk=startThisChunk+chunkSize

    print('   done in',time.time()-tic)

    tic=time.time()
    print('Starting to combine inBnotA')

    #now we need to find all (nxFrom,nyFrom) points that are in B but
    #not A, and the nxFrom/nyFrom/nxTo/nyTo/numTo lines for those
    #starting points to A. This is a straightforward append.
    indxInBnotA=array([ (p in inBnotA) for p in zip(nxFromB,nyFromB)])
    A['nxFrom'].append(B['nxFrom'].vindex[indxInBnotA],axis=0)                    
    A['nyFrom'].append(B['nyFrom'].vindex[indxInBnotA],axis=0)                    
    A['nxTo'].append(B['nxTo'].vindex[indxInBnotA],axis=0)                    
    A['nyTo'].append(B['nyTo'].vindex[indxInBnotA],axis=0)                    
    A['numTo'].append(B['numTo'].vindex[indxInBnotA],axis=0)                    

    print('   done in',time.time()-tic)
    

    return None


#so lets test this code
if __name__=="__main__":


    #==========================================================================================
    if False: #calculate connectivity matrix
        #turn of to avoid warnings, expriment for speed
        dask.config.set({'array.slicing.split_large_chunks': False})

        dataInFile='dataPaths/Faster_try1.zarr' #three months
        dataInFile=('../makeCommunityConnectivity/dataPaths/'+
                    'theAmericas_wholeGlobe_range100km_depthFrom_1m_to_500m_habitatTree_months01_to_12_starts_1m/2007.zarr')
        data=zarr.open(dataInFile)

        #what ages to get
        ageMin=13.0*86400 ; ageMax=13.0*86400
        ageMin=13.0*86400 ; ageMax=17.0*86400
        ages=data['age'][0,:] #assume all drifters have same age!
        ageIndex=logical_and(ages>=ageMin,ages<=ageMax)
        print('getting connectivity for PLD in days:',ages[ageIndex]/8.64e4)

        #filename to store connectivity data
        dirOut='dataPaths/jnkConnect.zarr'
        
        #now get a connectivty matrix for whole file
        if os.path.exists(dirOut):
            assert False,dirOut+' exists, delete by hand before running this'
        else:
            makeConnectFromFile(data,ageIndex,dirOut,minDay=-1,maxDay=30)

    #==========================================================================================
    if False: #calculate connectivity matrix which MIGHT only includes
             #points that start and land within a set of points, and
             #only include the shortest PLD that reaches the defined
             #habitat. See flag onlySettleInside to see if particles
             #that settle outside of habitat are included

        #include turn of to avoid warnings, expriment for speed
        dask.config.set({'array.slicing.split_large_chunks': False})

        #dataInFile='dataPaths/Faster_try1.zarr' #three months
        dataInFile='dataPaths/theAmericas_wholeGlobe_range100km_depth500m_habitatTree_months01_to_12_fixed_1m_2007.zarr' #full year
        dataInFile=('../makeCommunityConnectivity/dataPaths/'+
                    'theAmericas_wholeGlobe_range100km_depthFrom_1m_to_500m_habitatTree_months01_to_12_starts_1m/2007.zarr')
        data=zarr.open(dataInFile)

        #what ages to get
        ageMin=13.0*86400 ; ageMax=13.0*86400
        ageMin=13.0*86400 ; ageMax=14.0*86400
        ageMin=13.0*86400 ; ageMax=13.0*86400
        #ageMin=14.0*86400 ; ageMax=16.0*86400
        ages=data['age'][0,:] #assume all drifters have same age!
        ageIndex=logical_and(ages>=ageMin,ages<=ageMax)
        print('getting connectivity for PLD in days:',ages[ageIndex]/8.64e4)

        #make a dictionary of habitat where particles can start or
        #settle. This is just a set whose members are the points in
        #the habitat defined by (nx,ny) tuples. The code below creates
        #them from all starting points in the dataInFile. This is
        #horribly inefficient because many release points are
        #duplicated.
        if False:
            print('making habitat')
            firstTime=data['time'][0,0]
            indxFirstTime=data['time'][:,0]
            indxFirstTime=indxFirstTime==firstTime
            nxFrom=data['nx'][:,0]
            nyFrom=data['ny'][:,0]
            nxFrom=nxFrom[indxFirstTime]
            nyFrom=nyFrom[indxFirstTime]
            goodHabitat=set()
            goodHabitat.update(zip(nxFrom,nyFrom))
            print('   Done making good habitat')
            #save('jnk_goodHabitat.npz',goodHabitat)
        elif True:
            habitatFile='dataMatrices/theAmericas_wholeGlobe_range100km_depthFrom_1m_to_500m_habitatTree_months01_to_12_starts_1m/habitat_full.pckle'
            print('got goodHabitat from file ',habitatFile)
            goodHabitat=load(habitatFile,allow_pickle=True)
        else:
            print('use all habitat places within 3 grid points of land')
            print('making habitat')
            firstTime=data['time'][0,0]
            indxFirstTime=data['time'][:,0]
            indxFirstTime=indxFirstTime==firstTime
            nxFrom=data['nx'][:,0]
            nyFrom=data['ny'][:,0]
            nxFrom=nxFrom[indxFirstTime]
            nyFrom=nyFrom[indxFirstTime]
            allHabitat=set()
            allHabitat.update(zip(nxFrom,nyFrom))
            print('   Done making good habitat')
            nxVec=array([p[0] for p in allHabitat])
            nyVec=array([p[1] for p in allHabitat])
            gridDistDict=getGridDistanceDict(nxVec,nyVec)
            goodHabitat=set([k for k in gridDistDict if gridDistDict[k]<=3.0])
            
        #filename to store connectivity data
        dirOut='dataPaths/jnkConnect_limited.zip'

        #now get a connectivty matrix for whole file
        if os.path.exists(dirOut):
            assert False,dirOut+' exists, delete by hand before running this'
        else:
            makeConnectFromFileWithinHabitat(data,ageIndex,goodHabitat,dirOut,minDay=-1,maxDay=30.0,onlySettleInside=False)

        #===================================================================================
        #now plot and interact with it, if you wish

        #get connectivity
        dataConnect=zarr.open('dataPaths/jnkConnect_limited.zip','r')
        print('done loading points')

        #name of data file to get grid info from. %s will be replaced by the appropriate feild name.
        dataFileName='/data/break/pringle/mercatorFinal/2007/tmerc_phy_732716.nc'
        dataIn=xr.open_dataset(dataFileName)
        nav_lat=dataIn['nav_lat'][:]
        nav_lon=dataIn['nav_lon'][:]

        #get location of launch points
        nx=dataConnect['nxFrom'][:]
        ny=dataConnect['nyFrom'][:]
        latPnts=nav_lat[ny,nx]
        lonPnts=nav_lon[ny,nx]

        #click on point, see where children go.
        clf()
        #plot(nx,ny,'r.',zorder=2)
        #draw()

        if True:
            jnkNx=[p[0] for p in goodHabitat]
            jnkNy=[p[1] for p in goodHabitat]
            plot(jnkNx,jnkNy,'c,')
            #assert False,'think!'

        # title('cx for all goodHabitat points, r+, those that setelled')
        
        isFirst=True
        keepGoing=True #this is a hack to allow us to run just one point
        while True and keepGoing:

            #get starting point, and delete last plot if not first time
            print('click on point')
            if True:
                pntIn=ginput(n=1,timeout=-1)
            else:
                #for debugging, it is often useful to run same point over and over again...
                #empty points are nx= [2726 2726 2831] and ny= [2976 2977 1554]
                pntIn=[(2518.564438956629, 1932.966647608157)]
                keepGoing=False
            if not isFirst:
                han2.pop().remove() #for lines
                #han.pop().remove() #for lines
                han.remove() # for scatter
            isFirst=False

            #get locations of where to
            nxPnt=int(round(pntIn[0][0])); nyPnt=int(round(pntIn[0][1]))
            row=argmax(logical_and(nxPnt==nx,nyPnt==ny))
            nxPnts=dataConnect['nxTo'][row]
            nyPnts=dataConnect['nyTo'][row]
            numPnts=dataConnect['numTo'][row]
            try:
                print('  There are %d drifters from this point'%(sum(numPnts),))
            except:
                print('  Nothing found')
                
            #plot
            #han=plot(nxPnts,nyPnts,'r.')
            han=scatter(nxPnts,nyPnts,s=10,c=numPnts,zorder=5)
            han2=plot(nxPnt,nyPnt,'g*')

            #what region to zoom into
            #axis([2536,2834,2018,2266]) #SS GoM
            axis([2461,2749,1813,2148]) #East Coast

            draw()
            show()
            pause(0.1)

    #==========================================================================================
    elif False:
        #excersize routines for getting depth on model grid points

        #get connectivity
        data=zarr.open('dataPaths/jnkConnect.zarr','r')
        nxFrom=data['nxFrom'][:]
        nyFrom=data['nyFrom'][:]
        print('Making depth dict')
        depthDict=getDepthDict(nxFrom,nyFrom)
        print('   done')
        
        clf()
        scatter(nxFrom,nyFrom,s=5,c=[depthDict[(p[0],p[1])] for p in zip(nxFrom,nyFrom)],
                vmin=0.0,vmax=150)
        colorbar()

        draw()
        show()

        
    #==========================================================================================
    elif False:
        #excersize routines for getting distance from shore on model grid points

        #get connectivity
        data=zarr.open('dataPaths/jnkConnect.zarr','r')
        nxFrom=data['nxFrom'][:]
        nyFrom=data['nyFrom'][:]
        print('Making depth dict')
        distDict=getDistanceDict(nxFrom,nyFrom)
        print('   done')
        
        clf()
        scatter(nxFrom,nyFrom,s=5,c=[distDict[(p[0],p[1])] for p in zip(nxFrom,nyFrom)],vmin=0.0,vmax=100.0)
        colorbar()

        draw()
        show()

    #==========================================================================================
    elif False:
        #excersize routines for getting grid-distance from shore on model grid points

        #get connectivity
        data=zarr.open('dataPaths/jnkConnect.zarr','r')
        nxFrom=data['nxFrom'][:]
        nyFrom=data['nyFrom'][:]
        print('Making depth dict')
        gridDistDict=getGridDistanceDict(nxFrom,nyFrom)
        print('   done')
        
        clf()
        scatter(nxFrom,nyFrom,s=5,c=[gridDistDict[(p[0],p[1])] for p in zip(nxFrom,nyFrom)],vmin=0.0,vmax=5.0)
        colorbar()

        draw()
        show()
    #==========================================================================================
    elif False:
        #excersize code to trim connectivity matrix to only include
        #points that meet certain depth and distance criteria

        #get input connectivity data
        connectivityIn=zarr.open('dataPaths/jnkConnect.zarr','r')
        nxFrom=connectivityIn['nxFrom'][:]
        nyFrom=connectivityIn['nyFrom'][:]

        if False:
            print('Making depth dict')
            depthDict=getDepthDict(nxFrom,nyFrom)
            print('   done')
            depthMin=0.0; depthMax=100.0
            connectivityOut=trimConnectivity(connectivityIn,depthDict,depthMin,depthMax)
        else:
            print('Making dist dict')
            distDict=getGridDistanceDict(nxFrom,nyFrom)
            print('   done')
            distMin=0.0; distMax=2.1
            connectivityOut=trimConnectivity(connectivityIn,distDict,distMin,distMax)


        #if profiling, stop here
        #assert False,'done and think'
            
        #save it to disk
        import shutil
        try:
            shutil.rmtree('jnk.out')
        except:
            print('oops, jnk.out did not exist, you can ignore this')
        storeOut=zarr.DirectoryStore('jnk.out')
        rootOut=zarr.group(store=storeOut)
        zarr.copy_all(connectivityOut,rootOut)

        #read in what we wrote out, and run a simple model with it
        data=zarr.open('jnk.out')

        #name of data file to get grid info from. %s will be replaced by the appropriate feild name.
        dataFileName='/data/break/pringle/mercatorFinal/2007/tmerc_phy_732716.nc'
        dataIn=xr.open_dataset(dataFileName)
        nav_lat=dataIn['nav_lat'][:]
        nav_lon=dataIn['nav_lon'][:]

        #get location of launch points
        nx=data['nxFrom'][:]
        ny=data['nyFrom'][:]
        latPnts=nav_lat[ny,nx]
        lonPnts=nav_lon[ny,nx]

        #click on point, see where children go.
        clf()
        plot(nx,ny,'k,')
        draw()

        isFirst=True
        while True:

            #get starting point, and delete last plot if not first time
            print('click on point')
            pntIn=ginput(n=1,timeout=-1)
            if not isFirst:
                han2.pop().remove() #for lines
                #han.pop().remove() #for lines
                han.remove() # for scatter
            isFirst=False

            #get locations of where to
            nxPnt=int(round(pntIn[0][0])); nyPnt=int(round(pntIn[0][1]))
            row=argmax(logical_and(nxPnt==nx,nyPnt==ny))
            nxPnts=data['nxTo'][row]
            nyPnts=data['nyTo'][row]
            numPnts=data['numTo'][row]

            #plot
            #han=plot(nxPnts,nyPnts,'r.')
            han=scatter(nxPnts,nyPnts,s=5,c=numPnts)
            han2=plot(nxPnt,nyPnt,'g*')

            #what region to zoom into
            #axis([2536,2834,2018,2266]) #SS GoM
            axis([2461,2749,1813,2148]) #East Coast

            draw()
            show()
            pause(0.1)

    #==========================================================================================
    if False: #This code makes connectivity for several years. These
             #will then be used to test the append connectivity
             #matrices codes

        #include turn of to avoid warnings, expriment for speed
        dask.config.set({'array.slicing.split_large_chunks': False})

        #now make connectivity for several years
        for year in arange(2007,2015):
            #dataInFile='dataPaths/Faster_try1.zarr' #three months
            #dataInFile='../makeCommunityConnectivity/dataPaths/theAmericas_wholeGlobe_range100km_depth500m_habitatTree_months01_to_12_fixed_1m_%d.zarr'%(year,)
            dataInFile=('../makeCommunityConnectivity/dataPaths/'+
                        'theAmericas_wholeGlobe_range100km_depthFrom_1m_to_500m_habitatTree_months01_to_12_starts_1m/%d.zarr'%(year,))
            
            data=zarr.open(dataInFile)

            #what ages to get
            ageMin=14.0*86400 ; ageMax=14.0*86400
            ages=data['age'][0,:] #assume all drifters have same age!
            ageIndex=logical_and(ages>=ageMin,ages<=ageMax)
            print('getting connectivity for PLD in days:',ages[ageIndex]/8.64e4)

            #get a habitat file
            habitatFile='dataMatrices/theAmericas_wholeGlobe_range100km_depthFrom_1m_to_500m_habitatTree_months01_to_12_starts_1m/habitat_full.pckle'
            print('got goodHabitat from file ',habitatFile)
            goodHabitat=load(habitatFile,allow_pickle=True)

            #filename to store connectivity data
            dirOut='dataMatrices/jnkConnect_year_%d.zip'%(year,)

            #now get a connectivty matrix for whole file
            if os.path.exists(dirOut):
                shutil.rmtree(dirOut)

            tic=time.time()
            print('starting to make year',year)
            makeConnectFromFileWithinHabitat(data,ageIndex,goodHabitat,dirOut,
                                             minDay=-1,maxDay=30,onlySettleInside=False)
            print('done in',time.time()-tic)
            print(' ')

    #==========================================================================================
    if False: #test the code to create empty matrices and to combine matrices

        print('Making empty connectivity matrix')
        store=zarr.MemoryStore()
        A=makeEmptyConnectivity(store)
        #print(A['nxTo'].info)

        print('\ncombining one matrix to empty')
        tic=time.time()
        B=zarr.open('dataMatrices/jnkConnect_year_2007.zip','r')
        memBstore=zarr.MemoryStore() #copy into memory for speed
        memB=zarr.group(store=memBstore)
        zarr.copy(B,memB,name='/')
        combineConnectivity(A,memB)
        #print(A['nxTo'].info)
        print('In first thousand, the sum of numTo in A is',sum([sum(p) for p in A['numTo'][:]]))
        print('In first thousand, the sum of numTo in B is',sum([sum(p) for p in memB['numTo'][:]]))
        print('The total number of starting particles is',A.attrs['numberOfStartingTimes'])
        print('   done with combineConnectivity in',time.time()-tic)
        del memB
        del memBstore

        print('\ncombining second matrix to another with same starting locations')
        tic=time.time()
        B=zarr.open('dataMatrices/jnkConnect_year_2008.zip','r')
        memBstore=zarr.MemoryStore() #copy into memory for speed
        memB=zarr.group(store=memBstore)
        zarr.copy(B,memB,name='/')
        combineConnectivity(A,memB)
        #print(A['nxTo'].info)
        print('In first thousand, the sum of numTo is',sum([sum(p) for p in A['numTo'][:]]))
        print('In first thousand, the sum of numTo in B is',sum([sum(p) for p in memB['numTo'][:]]))
        print('The total number of starting particles is',A.attrs['numberOfStartingTimes'])
        print('   done with combineConnectivity in',time.time()-tic)
        del memB
        del memBstore

        print('\ncombining second matrix to another with same starting locations')
        tic=time.time()
        B=zarr.open('dataMatrices/jnkConnect_year_2009.zip','r')
        memBstore=zarr.MemoryStore() #copy into memory for speed
        memB=zarr.group(store=memBstore)
        zarr.copy(B,memB,name='/')
        combineConnectivity(A,memB)
        #print(A['nxTo'].info)
        print('In first thousand, the sum of numTo is',sum([sum(p) for p in A['numTo'][:]]))
        print('In first thousand, the sum of numTo in B is',sum([sum(p) for p in memB['numTo'][:]]))
        print('The total number of starting particles is',A.attrs['numberOfStartingTimes'])
        print('   done with combineConnectivity in',time.time()-tic)
        del memB
        del memBstore

        # print('\ncombining second matrix to another with same starting locations')
        # tic=time.time()
        # B=zarr.open('dataMatrices/jnkConnect_year_2010.zip','r')
        # memBstore=zarr.MemoryStore() #copy into memory for speed
        # memB=zarr.group(store=memBstore)
        # zarr.copy(B,memB,name='/')
        # combineConnectivity(A,memB)
        # #print(A['nxTo'].info)
        # print('In first thousand, the sum of numTo is',sum([sum(p) for p in A['numTo'][:]]))
        # print('In first thousand, the sum of numTo in B is',sum([sum(p) for p in memB['numTo'][:]]))
        # print('   done with combineConnectivity in',time.time()-tic)
        # del memB
        # del memBstore


        
        

        
