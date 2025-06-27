from pylab import *
from numpy import *
import zarr
from collections import defaultdict,Counter
import makeConnectivityModule as mcm
import popMatrixManipulate as pmm
import os
import shutil
import getEZfateFromOSN

#This code saves the connectivity matrix for a given month and larval
#depth behavior for the climatology described in
#https://github.com/JamiePringle/EZfate. This matrix defines where
#particles that start at a particular location end up. The file that
#containes this data will start with the name "E"

#it also calculates the transpose of this matrix, which defines the
#origin location of the any particles that reach a location. The file
#that contains this data will start with the name "Etranspose"

#IMPORTANT: in the function below, you define a "habitat", a region
#from which particles can start. YOU MUST DECIDE if you only want to
#keep pathways that start and end within the habitat, or you want to
#include pathways that start within the habitat and end outside. If
#you look for the origin of particles that end up outside of the
#habitat, this code will ONLY return starting points that were within
#the habitat. You must think if this is what you want, or if you want
#to expand your habitat. Note that you must both define a polygon to
#include the habitat, and make sure you loop over one or more of the
#five large sub-regions which contain the data if they are within your
#polygon ('AsiaPacific', 'EuropeAfricaMiddleEast', 'theAmericas',
#'Antarctica'). If you want to define a very large or global polygon,
#note that longitude goes from -180 to 180 degrees

#The connectivity data and their transpose are saved with a filename
#defined below in a directory called "connectivity" in this directory.
#This directory is created if it does not exist.
if not os.path.exists('connectivity'):
    os.makedirs('connectivity')

#Define where data we download will be saved. The data is cached here,
#so that we don't keep downloading data. However, the directory can
#get rather large. Any file in it can be deleted, as it will be
#re-downloaded if necessary. If you want to re-download the data (for
#example, to get an updated climatology file), then the file must be
#deleted from this directory. This directory is created if it does not
#exist
OSNdataDir='EZfateData/communityConnectivityMatrices/'

#make the function that will calculate the connectivity matrix and its
#transpose. It takes as it arguement the time the Lagrangian particle
#will drift in days and what year to get. Right now, this is limited to 2 to 60 days, by
#interval of 2 days. Please define which Lagrangian pathways you want
#in the function below.
def makeConnectivityAndTranspose(driftTimeDays,whatYear):


    #first, name your habitat -- this is your free choice
    habitatName='CmaenasHab'

    #do you want to keep paths that start in the habitat, but then
    #land outside?  if keepOutside=False, only include particles that
    #remain inside the domain
    keepOutside=True

    #what are the months in which dispersal happens?
    inMonths=[5]

    #what depth are the larvae. Can be 1,10,20 or 40 meters
    depth=1

    #how long do particles drift -- this is passed into the
    #function in the parameter driftTimeDays
    minPLD=driftTimeDays; maxPLD=minPLD

    #what vertical behavior? "fixed" means fixed to a depth
    #"starts" allows the larvae to drift passively in all 3
    #dimensions
    vertBehavior='fixed'

    #how far from land in 1/12th of a degree gridpoints should the habitat extend?
    offshoreGridExtent=2.1*sqrt(2.0)

    #Define the habitat that you will be examining. 
    #an easy way to define a polygon of lat lons is with
    #https://www.mapsdirections.info/en/draw-route-google-maps/
    #and save the KML file, and copy the lat/lon points from the KML file
    #make a polygon that encompasses your habitat. 
    regionPoly=[(-73.300781,34.537106),
                (-79.013672,37.660994),
                (-77.695313,41.32217),
                (-73.740234,43.782001),
                (-72.597656,46.990558),
                (-67.412109,50.912558),
                (-59.150391,52.438432),
                (-51.328125,53.12947),
                (-39.814453,50.690368),
                (-44.296875,42.044194),
                (-61.699219,35.257955)
                ]
        

    #define the name of the connectivity matrices, and where they will be stored, below
    #rootFileName defines the suffix of the filename. The connectivty forward in time
    #is in the file which starts with E, and the connectivity backwards in time is in
    #the matrices named Etranspose. 
    rootFileName='_%s_year%d_depth%d_minPLD%d_maxPLD%d_months%d_to_%d.zarr'%(habitatName,whatYear,
                                                                      depth,minPLD,maxPLD,
                                                                      amin(inMonths),amax(inMonths))
    transposeFileName='connectivity/Etranspose'+rootFileName
    connectFileName='connectivity/E'+rootFileName

    #now download the data and create the connectivity matrices and its inverse

    #load a number of matrices to play with, and combine them to get a seasonal matrix
    store=zarr.MemoryStore()
    E=mcm.makeEmptyConnectivity(store)

    #DEFINE IN LOOP BELOW WHICH REGION OR REGIONS YOU WISH TO INCLUDE
    #IN YOUR HABITAT -- there is no need to download the data for a
    #region if that region is not in the polygon defined by
    #regionPoly.
    for regionName in ['theAmericas']:#['AsiaPacific', 'EuropeAfricaMiddleEast', 'theAmericas']:
        for month in inMonths:
            print('   working on month',month,'in region',regionName,'for driftTimeDays',driftTimeDays,flush=True)
            matInFile=(OSNdataDir+
                       '%s/%dm/%s/'%(regionName,depth,vertBehavior)+
                       'year_%d_month%2.2d_minPLD%2.2d_maxPLD%2.2d.zip'%(whatYear,month,minPLD,maxPLD))

            #matInFile is the location of a file that define the
            #model grid and Lagrangian connectivity from the
            #EZfate project, as described in
            #https://github.com/JamiePringle/EZfate. The function
            #getEZfateFromOSN.getFileFromOSN() takes as an
            #argument a pathway to the EZfate data on the Open
            #Storage Network S3 bucket, and returns the path to a
            #local file that contains the same data. Details as to
            #where and how this done, and where the data is
            #stored, can be found in the getEZfateFromOSN module.
            matInFile=getEZfateFromOSN.getFileFromOSN(matInFile)
            matIn=zarr.open(matInFile,'r')

            #if a habitat is defined with distance from land
            #defined as gridcells from land, then we need to get
            #distance from land dictionary
            gridDict=mcm.getGridDistanceDict(matIn.nxFrom[:],matIn.nyFrom[:],landThresh=max(2.1,depth))

            #trim to all points within gridRadius of land
            #keep all true points, even if they fall outside of To, so we can
            #latter accurately count the number of points launched. 
            print('   trimming by distance from land')
            gridRadius=offshoreGridExtent
            matIn=mcm.trimConnectivity(matIn,gridDict,0.0,gridRadius,keepAllTo=True)

            mcm.combineConnectivity(E,matIn)
    print('DONE reading in',driftTimeDays,'which has',E.nxFrom.shape[0],'points',flush=True)
    print(' ',flush=True)

    #now, add the numLaunched variable to E. This is the total number
    #of To variables for each from variable
    numLaunched=[sum(p) for p in E.numTo]
    root=zarr.group(store=zarr.MemoryStore())
    numLaunchedZarrVar=root.empty(shape=(len(numLaunched),),name='numLaunched',dtype='i4')
    numLaunchedZarrVar[:]=array(numLaunched)
    zarr.convenience.copy(numLaunchedZarrVar,E,name='numLaunched')

    #now trim to only include habitat inside of the regionPoly
    #polygon
    print('trimming by polygon')
    E=mcm.trimConnectivity_byPoly(E,regionPoly,keepInPoly=True,keepAllTo=keepOutside)

    if not keepOutside:
        #now trim connectivity again, but this time get rid of *To
        #points that land outside of habitat
        print('trimming To points')
        E=mcm.trimConnectivity(E,gridDict,0.0,gridRadius,keepAllTo=False)

    #write out connectivity matrix
    #breakpoint()
    print('writing connectivity to disk',driftTimeDays,flush=True)
    if os.path.exists(connectFileName):
        shutil.rmtree(connectFileName)
    EoutStore=zarr.DirectoryStore(connectFileName)
    Eout=zarr.group(store=EoutStore)
    zarr.copy(E,Eout,name='/')
    print('   done writing',driftTimeDays,flush=True)

    #now invert matrix so we have the connectivity from where the
    #larvae end up to where they started (the inverse of what we
    #just calculated)
    EtransposeStore=zarr.MemoryStore() #in memory
    Etranspose=pmm.invertConMatrix(E,EtransposeStore)

    # normalizes the numTo array, which is the number of
    # lagrangian particles which go to a given point in E by the
    # sum of numTo, to get the likelyhood that the point in
    # (nxFrom,nyFrom) went to one the (nxTo,nyTo) points in the
    # matrix.  This normalized array is put in the variable
    # propTo.  Also calculates the cumilative sum of this matrix
    # and places it in cumSumPropTo

    print('normalizing driftTimeDays',driftTimeDays,flush=True)
    pmm.normalizeInvertConMatrix(Etranspose)
    print('   done normalizing',driftTimeDays,flush=True)

    #now write out to disk
    print('writing transpose to disk',driftTimeDays,flush=True)
    if os.path.exists(transposeFileName):
        shutil.rmtree(transposeFileName)
    EoutStore=zarr.DirectoryStore(transposeFileName)
    Eout=zarr.group(store=EoutStore)
    zarr.copy(Etranspose,Eout,name='/')
    print('   done writing',driftTimeDays,flush=True)
    print('\n\n')

    return connectFileName,transposeFileName #return the fileNames of the matrices we made


#create connectivity matrices for a single drift duration.
#look in the plotting code to learn how to interpret this data
if __name__ == "__main__":
    #now run the function to make the habitat for a given lagrangian
    #drift duration (here called PLD) and several years
    whatPLD=40
    for whatYear in [2020,2021,2022]:
        connectFileName,transposeFileName=makeConnectivityAndTranspose(whatPLD,whatYear)

    #now load the connectivity matrices for the last year so you can look at them. 
    E=zarr.open(connectFileName,'r')
    Etranspose=zarr.open(transposeFileName,'r')
