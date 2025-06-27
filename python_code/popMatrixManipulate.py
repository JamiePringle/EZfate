from pylab import *
from numpy import *
import zarr
import numcodecs
from collections import defaultdict,Counter
import makeConnectivityModule as mcm
import time

#this module contains routines to manipulate the connectivity matrices
#from 01_makeConnectivity_byMonth_andYear.py and
#02_makeClimatological_matrices_OneMonth.py and
#02_runYearInParallel.com

#define the type of a model coordinate; must match that in makeConnectivityModule.py
coordType=uint16

# the type of the number of drifters that make a particular
# trip. Warning, danger of overflow.... Must match that in makeConnectivityModule.py
typeOfSum=uint16

#@profile
def invertConMatrix(E,EtransposeStore):
    '''invertConMatrix(E,EtransposeStore):

    Takes as input a zarr dataset E and the store for a zarr dataset
    EtransposeStore and returns the zarr dataset Etranspose (stored in
    EtransposeStore) with the following data:

    make an matrix transpose of a population connectivity matrix, so that
    nxTo is a grid index of type coordType, where the larvae ended up
    nyTo is a grid index of type coordType, where the larvae ended up
    nxFrom: grid indices of starting points of drifters, stored as ragged array
    nyFrom: grid indices of starting points of drifters, stored as ragged array
    numFrom: number of points that started in (nxFrom,nyFrom) that ended up in (nxTo,nyTo), stored as a ragged array

    returns Etranspose with data added

    '''

    # #figure out how many unique (nxTo,nyTo) pairs there are and keep
    # #track of the number of pairs. This can be surprizingly expensive,
    # print('finding unique to points')
    # allToPnts=set()
    # for nRow in range(len(E['nxFrom'])):
    #     allToPnts.update(zip(E['nxTo'][nRow],E['nyTo'][nRow]))
    # nTotal=len(allToPnts)
    # print('   done for %d unique(nxTo,nyTo) pairs'%(nTotal,))

    #THIS MIGHT BE ABSURDLY SLOW
    #ok, now collect the data. this will be done in memory. make a
    #dictionary whose key is (nxTo,nyTo) and whose value is a
    #counter. the value counter will have a key of
    #(nxFrom,nyFrom) and a value of the number of points that went
    #from (nxFrom,nyFrom) to (nxTo,nyTo)
    print('Making transform for ',len(E['nxFrom']),'rows'); tic=time.time()
    transDict=defaultdict(Counter)
    tic=time.time()

    #for efficiency, I found it very, very necessary to read nxTo,nyTo and numTo into memory.
    #reading it in a row at a time was much, much too slow
    this_nxTo=E['nxTo'][:]
    this_nyTo=E['nyTo'][:]
    this_numTo=E['numTo'][:]
    this_nxFrom=E['nxFrom'][:]
    this_nyFrom=E['nyFrom'][:]
    
    for nRow in range(len(E['nxFrom'])):
        if remainder(nRow,10000)==0:
            print('row',nRow,'or',nRow/len(E['nxFrom'])*100,'percent done,',
                  len(transDict),'unique (nxTo,nyTo) pairs','in',time.time()-tic,'s',flush=True)
            tic=time.time()
        for toPntAndCount in zip(this_nxTo[nRow],this_nyTo[nRow],this_numTo[nRow]):
            toPnt=toPntAndCount[:2]
            numToPnt=toPntAndCount[2]
            transDict[toPnt].update({(this_nxFrom[nRow],this_nyFrom[nRow]):
                                     numToPnt}) #update the count for tuple(nxFrom,nyFrom)
    nTotal=len(transDict)
    print('   done for %d unique (nxTo,nyTo) pairs'%(nTotal,),'in a time',time.time()-tic)

    #breakpoint()
    
    #create zarr dataset
    chunkSize=1000 #should experiment with this
    Etranspose=zarr.group(store=EtransposeStore)
    
    nxTo=Etranspose.empty(shape=(nTotal,),name='nxTo',dtype=coordType,chunks=chunkSize)
    nxTo[:]=array([p[0] for p in transDict])
    
    nyTo=Etranspose.empty(shape=(nTotal,),name='nyTo',dtype=coordType,chunks=chunkSize)
    nyTo[:]=array([p[1] for p in transDict])

    #create nxFrom, nyFrom and numFrom in zarr array
    nxFrom=Etranspose.empty(shape=(nTotal,),name='nxFrom',dtype=object,
                    chunks=chunkSize,object_codec=numcodecs.VLenArray(coordType))
    nyFrom=Etranspose.empty(shape=(nTotal,),name='nyFrom',dtype=object,
                    chunks=chunkSize,object_codec=numcodecs.VLenArray(coordType))
    numFrom=Etranspose.empty(shape=(nTotal,),name='numFrom',dtype=object,
                     chunks=chunkSize,object_codec=numcodecs.VLenArray(typeOfSum))

    #write nxFrom, nyFrom, numFrom
    tic=time.time()
    print('writting out nxFrom, nyFrom and numFrom')
    nRow=-1
    for pnt in transDict:
        nRow+=1 #leaning heavily hear on fact the dicts are ordered in python version>3.6

        jnkArray=zeros((1,),dtype=object) #an array of object of shape (1,)

        jnk=array([p[0] for p in transDict[pnt].keys()])
        jnkArray[0]=jnk
        nxFrom[nRow]=jnk

        jnk=array([p[1] for p in transDict[pnt].keys()])
        jnkArray[0]=jnk
        nyFrom[nRow]=jnk

        jnk=array([p for p in transDict[pnt].values()])
        jnkArray[0]=jnk
        numFrom[nRow]=jnk

    print('   done writing in',time.time()-tic)

    return Etranspose

def normalizeInvertConMatrix(Etranspose):
    '''normalizeInvertConMatrix(Etranspose):

    This routine normalizes the numFrom array in Etranspose by the sum
    of numFrom, which as described in in "BackwardsProp.pdf" makes
    numFrom the likelyhood that a point in (nxTo,nyTo) came from one
    the (nxFrom,nyFrom) points in the matrix.  This normalized array
    is put in the variable propFrom

    it also calculates the cumilative sum of this matrix and places it in cumSumPropFrom

    This routine takes Etranspose as an input, and returns no output
    since Etranspose is modified to include a new variable, propFrom.

    '''

    propFrom=Etranspose.empty(shape=Etranspose['numFrom'].shape,name='propFrom',dtype=object,
                              chunks=Etranspose['numFrom'].chunks,
                              object_codec=numcodecs.VLenArray(type(1.0)))

    cumSumPropFrom=Etranspose.empty(shape=Etranspose['numFrom'].shape,name='cumSumPropFrom',dtype=object,
                              chunks=Etranspose['numFrom'].chunks,
                              object_codec=numcodecs.VLenArray(type(1.0)))

    if False: #old, slow low memory way
        for nRow in range(len(Etranspose['numFrom'])):        
            propFrom[nRow]=Etranspose['numFrom'][nRow]/sum(Etranspose['numFrom'][nRow])
            cumSumPropFrom[nRow]=cumsum(Etranspose['numFrom'][nRow]/sum(Etranspose['numFrom'][nRow]))
    else:
        #faster, preallocate arrays, write at once
        propFromTemp=empty(Etranspose['numFrom'].shape,dtype=object)
        cumSumPropFromTemp=empty(Etranspose['numFrom'].shape,dtype=object)
        numFrom=Etranspose['numFrom'][:]
        tic=time.time()
        for nRow in range(len(numFrom)):
            if remainder(nRow,10000)==0:
                print('row',nRow,'or',nRow/len(numFrom)*100,'percent done in',time.time()-tic,'s',flush=True)
                tic=time.time()
            propFromTemp[nRow]=numFrom[nRow]/sum(numFrom[nRow])
            cumSumPropFromTemp[nRow]=cumsum(numFrom[nRow])/sum(numFrom[nRow])
        propFrom[:]=propFromTemp
        cumSumPropFrom[:]=cumSumPropFromTemp

    return None

def normalizeConMatrix(E):
    '''normalizeConMatrix(E):

    This routine normalizes the numTo array in E by the sum of numTo,
    which makes numTo the likelyhood that the point in (nxFrom,nyFrom)
    went to one the (nxTo,nyTo) points in the matrix.  This
    normalized array is put in the variable propTo

    it also calculates the cumilative sum of this matrix and places it in cumSumPropTo

    This routine takes E as an input, and returns no output
    since E is modified to include a new variable, propTo.

    '''

    propTo=E.empty(shape=E['numTo'].shape,name='propTo',dtype=object,
                              chunks=E['numTo'].chunks,
                              object_codec=numcodecs.VLenArray(type(1.0)))

    cumSumPropTo=E.empty(shape=E['numTo'].shape,name='cumSumPropTo',dtype=object,
                              chunks=E['numTo'].chunks,
                              object_codec=numcodecs.VLenArray(type(1.0)))

    if False: #old, slow, low memory
        for nRow in range(len(E['numTo'])):
            #print(nRow,'of',len(E['numTo']))
            propTo[nRow]=E['numTo'][nRow]/sum(E['numTo'][nRow])
            cumSumPropTo[nRow]=cumsum(E['numTo'][nRow]/sum(E['numTo'][nRow]))
    else:
        #faster, preallocate arrays, write at once
        propToTemp=empty(E['numTo'].shape,dtype=object)
        cumSumPropToTemp=empty(E['numTo'].shape,dtype=object)
        numTo=E['numTo'][:]
        for nRow in range(len(numTo)):
            #print(nRow,'of',len(numTo))
            propToTemp[nRow]=numTo[nRow]/sum(numTo[nRow])
            cumSumPropToTemp[nRow]=cumsum(numTo[nRow]/sum(numTo[nRow]))
        propTo[:]=propToTemp
        cumSumPropTo[:]=cumSumPropToTemp

    return None

    

if __name__=="__main__":
    import os

    #the function getEZfateFromOSN.getFileFromOSN() takes as an argument a pathway to the
    #EZfate data on the open storage network S3 bucket, and returns the path to a local
    #file that contains the same data. Details as to where and how this done, and where
    #the data is stored, can be found in the getEZfateFromOSN module.
    import getEZfateFromOSN


    inParam=30 #PLD in days    
    
    #what habitat, depth and months
    inMonths=arange(4,6) #two months of data
    depth=1
    minPLD=inParam; maxPLD=minPLD
    vertBehavior='starts'

    #habitatName is the name of the habitat that defines the Transposed matrix.
    habitatName='Global_2gp_test'

    #ok, if True, make a transposed matrix from data (and save a copy to
    #disk), if False, read it from disk...
    rootFileName='_%s_depth%d_minPLD%d_maxPLD%d_months%d_to_%d.zarr'%(habitatName,
                                                        depth,minPLD,
                                                        maxPLD,amin(inMonths),amax(inMonths))
    transposeFileName='transposes/Etranspose_test'+rootFileName
    connectFileName='transposes/E_test'+rootFileName

    if True: #make connectivity; if false, read existing file from disk
        #load a number of matrices to play with, and combine them to get a seasonal matrix
        store=zarr.MemoryStore()
        E=mcm.makeEmptyConnectivity(store)
        for regionName in ['AsiaPacific', 'EuropeAfricaMiddleEast', 'theAmericas']:
            for month in inMonths:
                print('working on month',month,'in region',regionName,'for inParam',inParam,flush=True)
                matInFile=('EZfateData/communityConnectivityMatrices/'+
                           '%s/%dm/%s/'%(regionName,depth,vertBehavior)+
                           'climatology_month%2.2d_minPLD%2.2d_maxPLD%2.2d.zip'%(month,minPLD,maxPLD))
                matInFile=getEZfateFromOSN.getFileFromOSN(matInFile)
                matIn=zarr.open(matInFile,'r')

                #get distance from land dictionary
                gridDict=mcm.getGridDistanceDict(matIn.nxFrom[:],matIn.nyFrom[:],landThresh=max(2.1,depth))

                #trim to all points within gridRadius of land
                print('trimming by distance from land')
                gridRadius=2.1*sqrt(2)
                matIn=mcm.trimConnectivity(matIn,gridDict,0.0,gridRadius)

                #trim to only include points in regionPoly
                #print('trimming by polygon')
                #matIn=mcm.trimConnectivity_byPoly(matIn,regionPoly,keepInPoly=True)
                #breakpoint()

                mcm.combineConnectivity(E,matIn)
        print('DONE reading in',inParam,'which has',E.nxFrom.shape[0],'points',flush=True)
        print(' ',flush=True)

        #write out connectivity matrix
        print('writing connectivity to disk',inParam,flush=True)
        if os.path.exists(connectFileName):
            shutil.rmtree(connectFileName)
        EoutStore=zarr.DirectoryStore(connectFileName)
        Eout=zarr.group(store=EoutStore)
        zarr.copy(E,Eout,name='/')
        print('   done writing',inParam,flush=True)
    else:
        #read in and copy to memory
        print('reading in E and putting it into memory')
        Ein=zarr.open(connectFileName,'r')
        store=zarr.MemoryStore()
        E=zarr.group(store=store)
        zarr.copy(Ein,E,name='/')
        print('   done')
                      


    print('starting to transpose')
    tic=time.time()
    #now invert matrix
    EtransposeStore=zarr.MemoryStore() #in memory
    Etranspose=invertConMatrix(E,EtransposeStore)
    print('    done in',time.time()-tic)

    #can I normalize it?
    tic=time.time()
    print('normalizing inParam',inParam,flush=True)
    normalizeInvertConMatrix(Etranspose)
    print('   done normalizing',inParam,'in',time.time()-tic,flush=True)

    #now write out to disk
    if False:
        print('writing transpose to disk',inParam,flush=True)
        if os.path.exists(transposeFileName):
            shutil.rmtree(transposeFileName)
        EoutStore=zarr.DirectoryStore(transposeFileName)
        Eout=zarr.group(store=EoutStore)
        zarr.copy(Etranspose,Eout,name='/')
        print('   done writing',inParam,flush=True)


