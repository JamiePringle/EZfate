
library(ncdf4)
library(tidyverse,quietly=TRUE,warn.conflicts=FALSE) #get it to shut up!
library(sf,quietly=TRUE,warn.conflicts=FALSE)
library(collections,quietly=TRUE,warn.conflicts=FALSE)

#if true, enable parallel operations. In general, you want to do this!
#it makes, among other things, the combination of connectivity data 
#much faster.
if (TRUE) {
  runParallel=TRUE
  library(furrr,quietly=TRUE,warn.conflicts=FALSE)
  future::plan(multisession)
  nProcessors=parallel::detectCores()
}else{
  runParallel=FALSE
}

#setup debugging output if TRUE
if (TRUE) {
  #setup for debugging
  debugOutput=TRUE
  library(tictoc,quietly=TRUE,warn.conflicts=FALSE)
} else {
  debugOutput=FALSE
}

#===============================================================================
#define routines to download data, and where the data is obtained from and where
#it will be stored. 
# This code downloads connectivity data IF IT HAS NOT ALREADY BEEN DOWNLOADED.
# If you wish to download files which have been updated on the server, you 
# must delete them from the local directory first. 
# 
# To download the data you need to define several parameters:
#
# regionName: a string that defines what region to download. For example,
# "theAmericas", "theArctic", etc
#
# year: a number that specifies the year of the particle tracks OR you must
# specify that you want the climatological average over all years.
#
# depth: a number that specifies the depth the particles are released at in
# meters. There are a finite number of depths (e.g. 1, 20 and 40).
#
# verticalBehavior: a string that specifies the vertical behavior of the
# particles as they are moved by the modeled currents. E.g. "fixed" or "starts".
#
# minPLD and maxPLD: an integer that specifies the beginning and end of the time
# the particle can enter suitable habitat. Often will be the same number. 
# 
# month: and integer which specifies the month of the data

# We need to specify the root URL which contains the data
rootDataURL<-'http://oxbow.sr.unh.edu/data/RcommunityConnectivityMatrices/'

# and specify the layout of the files below rootDataURL; this string will be
# processed by sprintf, and as currently defined expects
# (depth,verticalBehavior, a string that is either year_2007 (or appropriate
# year) or "climatological", month, minPLD, maxPLD). Use file.path so that the
# path is appropriately made on all architectures.  
dataLayout<-file.path('%s','allPoints','%dm','%s','%s_month%2.2d_minPLD%2.2d_maxPLD%2.2d.RDS')

getConnectivityData<-function(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD,dataDir='connectivityData'){
  #this code downloads into dataDir the connectivity data specified by the arguements
  #regionName,depth,year,verticalBehavior,month,minPLD and maxPLD, as defined in the comments 
  #above. 
  #
  #Note that year can either be an integer (e.g. 2007) or the string "climatological"
  #depending on what is desired. 
  #
  #NOTE WELL -- DATA WILL NOT BE RE-DOWNLOADED IF THE FILE ALREADY EXISTS
  #
  #the function returns the data.frame that has been downloaded. 
  
  #first, make a string to encode year information or climatology, as needed
  if (year=='climatology'){
    yearStr<-'climatology'
  } else {
    yearStr<-sprintf('year_%d',year)
  }
  
  #put download path and filename together, along with URL of data
  downloadPath<-sprintf(dataLayout,regionName,depth,verticalBehavior,yearStr,month,minPLD,maxPLD)
  downloadToFile<-file.path(dataDir,downloadPath)
  dataURL<-paste(rootDataURL,downloadPath,sep='')
  
  #only make True for debugging
  if (FALSE){
    #print out debuging information
    print(dataLayout)
    print(downloadToFile)
    print(basename(downloadToFile))
    print(dirname(downloadToFile))
    print(dataURL)
  }
  
  #make directory for data
  dirName<-dirname(downloadToFile)
  if (!dir.exists(dirName)){
    print(paste('creating directory',dirName))
    dir.create(dirName,recursive=TRUE)
  }
  
  #download file to downloadToFile if it does not exist
  if (file.exists(downloadToFile)) {
    print(paste('file',downloadToFile,'exists, no download'))
  } else{
    options(timeout = max(3000, getOption("timeout")))
    download.file(dataURL,downloadToFile,mode='wb')
  }
  
  #load the data.frame and return it
  E<-readRDS(downloadToFile)
  return(E)
}

getGridData<-function(dataDir='connectivityData'){
  #this code gets the model grid data file. It has one arguement, dataDir, 
  #which should be the same as the argument to getConnectivityData(). It returns
  #the variable associated with the netCDF file. 
  dirName<-file.path(dataDir,'EZfateFiles')
  downloadToFile<-file.path(dataDir,'EZfateFiles','model_depth_and_distance.nc')
  dataURL<-paste(rootDataURL,'EZfateFiles/model_depth_and_distance.nc',sep='')
  
  if (!dir.exists(dirName)){
    print(paste('creating directory',dirName))
    dir.create(dirName,recursive=TRUE)
  }
  
  #download file to downloadToFile if it does not exist
  if (file.exists(downloadToFile)) {
    print(paste('file',downloadToFile,'exists, no download'))
  } else{
    options(timeout = max(3000, getOption("timeout")))
    download.file(dataURL,downloadToFile,mode='wb')
  }
  
  data<-nc_open(downloadToFile)
  return(data)
}


#===============================================================================
#now start by loading into memory the model grid data. this enables functions that
#takes vectors of nx,ny and returns vectors of lon,lat from the Mercator GLORYS
#1/12 degree run rho grid. It also enables the connectivity data to be subset by
#water depth or distance from land.

#use the getGridData() function to download and access the grid data.
data<-getGridData()

lon<-ncvar_get(data,'nav_lon')
lat<-ncvar_get(data,'nav_lat')
gridDepth<-ncvar_get(data,'depth') #the ocean depth on model grid
gridDist<-ncvar_get(data,'gridDist') #the distance from land in gridCells
gridDistKm<-ncvar_get(data,'dist') #the distance from land in kilometers

#shut dplyr.summarise up...
options(dplyr.summarise.inform = FALSE) #when using summarize below, it gets _way_ to chatty

#===============================================================================
#the following functions are utility functions to manipulate the connectivity
#data.frames

#this takes in vectors of nx and ny and returs a list with lonVec, latVec,
#depthVec, the latitude, longitude and depth of the nx,ny points. BECAUSE THE
#INDICES COME FROM PYTHON/C I NEED TO ADD 1 TO THEM? CODE OFFSET TO INDICES AS
#NOFF
Noff=1
nxny2lonlat <- function(nx,ny) {
  latVec<-lat[cbind(nx+Noff,ny+Noff)]
  lonVec<-lon[cbind(nx+Noff,ny+Noff)]
  depthVec<-gridDepth[cbind(nx+Noff,ny+Noff)]
  returnVal<-list("latVec"=latVec, "lonVec"=lonVec,'depthVec'=depthVec)
  return(returnVal)
}

#just return lat, or just return lon, for use in following function
nxny2lon <- function(nx,ny) {
  lonVec<-lon[cbind(nx+Noff,ny+Noff)]
  return(lonVec)
}
nxny2lat <- function(nx,ny) {
  latVec<-lat[cbind(nx+Noff,ny+Noff)]
  return(latVec)
}
#this function takes a dispersal structure and returns its with all nx* and ny* variables 
#suplemented with lon* and lat* values
addLatLon<-function(E){
  returnVals<-nxny2lonlat(E$nxFrom,E$nyFrom)
  
  #do the "from", each a simple vector
  E$lonFrom<-returnVals$lonVec
  E$latFrom<-returnVals$latVec
  
  # #apply to each row for lonTo and latTo
  # print('start apply 1')
  # E$lonTo<-apply(E,1,nxny2lon)
  # print('end apply 2')
  # E$latTo<-apply(E,1,nxny2lat)
  # print('done with applies')

  # this is not memory efficient, but is it fast?
  nxToList<-apply(E,1,function(a) {return(a$nxTo)})
  nyToList<-apply(E,1,function(a) {return(a$nyTo)})
  E$lonTo<-map2(nxToList,nyToList,nxny2lon)
  E$latTo<-map2(nxToList,nyToList,nxny2lat)
  
  return(E)
}

#===============================================================================
#the following functions subset the connectivity data into smaller regions

trimToOnlyFromPlaces<-function(EplusTrimmed){
  #this code trims numTo,nxTo,nyTo,lonTo,latTo so they only contain (nx,ny)
  #points that are (nxFrom,nyFrom) points
  #make dictionary of (nxFrom,nyFrom), so we can quickly check if point is in
  #the starting points
  startPoints<-dict(EplusTrimmed$nxFrom,map2(EplusTrimmed$nxFrom,EplusTrimmed$nyFrom,c)) #items, keys!
  
  #now, lets see if for this start point, which nxTo,nyTo are in the startPoints
  #then trim nxTo,nyTo and numTo to only include those points
  for (np in 1:nrow(EplusTrimmed)){
    #first find which *To's are also starting points
    nxTo<-EplusTrimmed$nxTo[[np]]
    nyTo<-EplusTrimmed$nyTo[[np]]
    numTo<-EplusTrimmed$numTo[[np]]
    allTo<-map2(nxTo,nyTo,c)
    inStarts<-unlist(map(allTo,startPoints$has)) #map returns list of lists, unlist makes into a single list
    
    #print(np)
    #now put trimmed records back into EplusTrimmed
    EplusTrimmed$nxTo[[np]]<-nxTo[inStarts]
    EplusTrimmed$nyTo[[np]]<-nyTo[inStarts]
    EplusTrimmed$numTo[[np]]<-numTo[inStarts]
    EplusTrimmed$lonTo[[np]]<-EplusTrimmed$lonTo[[np]][inStarts]
    EplusTrimmed$latTo[[np]]<-EplusTrimmed$latTo[[np]][inStarts]
  }
  return(EplusTrimmed)
}

subsetConnectivity_byPolygon<-function(Eplus,limitPoly,trimTo=FALSE){
  #this code takes a connectivity dataframe that has had lat and lon added to it
  #with addLatLon, and then subsets it to only include points that start within
  #limitPoly. If trimTo=TRUE, it will not include To points that are not at
  #locations which have nxFrom,nyFrom points. If trimTo=FALSE (the default), it
  #will include all To points.
  
  #make an sf sfc object (see https://geocompr.robinlovelace.net/index.html) of all the *From points
  thePoints<-st_multipoint(matrix(cbind(Eplus$lonFrom,Eplus$latFrom),nrow(Eplus),2))
  thePoints=st_sfc(thePoints,crs=4326)
  
  #make into list of points
  thePoints<-st_cast(thePoints,'POINT')
  
  #which points are in limitPoly
  inLimit<-st_contains(limitPoly,thePoints)
  
  #only include nxFrom,nyFrom that are inside polygon
  EplusTrimmed<-Eplus[inLimit[[1]],]
  
  #if trimTo=True, then get all possible starting points from EplusTrimmed and
  #store starting points as (nx,ny) in a dictionary from the collections
  #library. use this dictionary to quickly determine find the
  #nxTo,nyTo,lonTo,latTo and numTo which fall within starting region, and only
  #include those in the output.
  if (trimTo) {
    EplusTrimmed<-trimToOnlyFromPlaces(EplusTrimmed)
  }
  
  return(EplusTrimmed)
}

subsetConnectivity_byGridValue<-function(Eplus,values,minValue,maxValue,trimTo=FALSE){
  #this code takes a connectivity dataframe and a matrix of values on the model
  #grid and subsets the connectivity to only include points whose value in
  #values lie withing minValue<=values[Eplus$nxFrom,Eplus$nyFrom]<maxValue. If
  #trimTo=TRUE, it will not include To points that are not at locations which
  #have nxFrom,nyFrom points. If trimTo=FALSE (the default), it will include all
  #To points.
  
  #make a vector of values[Eplus$nxFrom,Eplus$nyFrom]. NOTE WELL, since nx and
  #ny in Eplus were computed in python, we must add 1 to the values before
  #looking them up in R
  theValues<-values[cbind(Eplus$nxFrom+1,Eplus$nyFrom+1)]
  
  #index true for points we want to keep
  indx<-(theValues>=minValue)  & (theValues<maxValue)
  
  #only include nxFrom,nyFrom that are inside polygon
  EplusTrimmed<-Eplus[indx,]
  
  #if trimTo=True, then get all possible starting points from EplusTrimmed and
  #store starting points as (nx,ny) in a dictionary from the collections
  #library. use this dictionary to quickly determine find the
  #nxTo,nyTo,lonTo,latTo and numTo which fall within starting region, and only
  #include those in the output.
  if (trimTo) {
    EplusTrimmed<-trimToOnlyFromPlaces(EplusTrimmed)
  }
  
  return(EplusTrimmed)
}

#===============================================================================
#the following functions combine two connectivity data.frames into a single data
#frame
mergeColumnByRowGroupsLemma<-function(x){
  #this code is designed to be called by summarize() on a column in a few rows
  #that have been grouped by group_by(), and combine the list in the columns. This
  #merges the lists in nxTo,nyTo and numTo. 
  xout<-list(c(x,recursive=TRUE))
  #print('stop here')
  return(xout)
}

deDuplicateToDataLemma<-function(nxTo,nyTo,numTo) {
  #after merging the two connectivity datas in combineConnectivityData(), there
  #are multiple identicle (nxTo,nyTo) pairs, along with the associated numTo's
  #
  #this code takes as an argument a row from the combined data.frame, called
  #aRow, quickly finds all common (nxTo,nyTo) pairs in that row, adds their
  #associated numTo's together, and returns the updated row as aRowOut
  
  #this code is a bit tricky. Use a dictionary d whose keys will be c(nxTo,nyTo)
  #and whose value is numTo. Use map to populate dictionary by adding numTo
  #for each (nxTo,nyTo) pair to the entry for the dictionary, taking advantage
  #of the fact that you can set a default value for the dictionary if they 
  #key is not in the dictionary. The output of the pmap is not used; the 
  #whole point is the population of the dictionary. 
  
  if (TRUE) {
    #first Try, with dict 2.9 ms/row
    d<-ordered_dict()
    pmap(list(nxTo,nyTo,numTo),function(a,b,c) (d$set(c(a,b),d$get(c(a,b),0)+c)))
    nxToOut<-c(map(d$keys(),function(a)(a[1])),recursive=TRUE)
    nyToOut<-c(map(d$keys(),function(a)(a[2])),recursive=TRUE)
    numToOut<-c(map(d$values(),function(a)(a)),recursive=TRUE)
  }
  if (FALSE){
    #data.frame, group_by,summarize 2.7 ms/row
    #but is not right? Check with sum(E[1,]$numTo[[1]])
    outData<- data.frame(nxTo=nxTo,nyTo=nyTo,numTo<-numTo) %>%
      group_by(nxTo,nyTo) %>% summarise(numTo=sum(numTo))
    nxToOut<-outData$nxTo
    nyToOut<-outData$nyTo
    numToOut<-outData$numTo
  }

  
  #print('hi')
  return(list(nxToOut,nyToOut,numToOut))
}

combineConnectivityData<-function(E1,E2){
  #this code takes two connectivity data.frames and combines them. Note -- if
  #they have had lat and lon added to them, it will be lost. So you will need to
  #add lat and lon to them again.
  
  #the process below removes the lat,lon variables. So lets check if they 
  #exist so we can add them later
  haveLatLon <- ("lonFrom" %in% colnames(E1)) | ("lonFrom" %in% colnames(E2))
  
  #now combine by common nxFrom,nyFrom, and concatenate nxTo, nyTo and numTo
  #
  #however, this has an issue. In each row of the data frame, for nxTo, nyTo and
  #numTo there are data for common (nxTo,nyTo) pairs. Note -- not all sets of
  #(nxTo,nyTo) data exist in pairs. To combine these, use a function defined
  #above called combineToDataLemma(nxTo,nyTo,numTo) which takes in the data for
  #nxTo,nyTo and numTo for each row, and combines them. Use map to iterate this
  #function over all the rows in the data  
  jnk<-rbind(E1,E2) #make a single frame with all data
  combinedE<-jnk %>% group_by(nxFrom,nyFrom) %>% summarise(nxTo=mergeColumnByRowGroupsLemma(nxTo),
                                                      nyTo=mergeColumnByRowGroupsLemma(nyTo),
                                                      numTo=mergeColumnByRowGroupsLemma(numTo),numLaunched=sum(numLaunched)) 
  
  #now this will go through each of nxTo,nyTo and numTo entries, and combine
  #common (nxTo,nyTo) pairs, along with the sum of their numTo's and return
  #it as a list of lists of the (nxTo,nyTo,numTo) data
  if (debugOutput){tic(msg='Time to combine connectivity data:')}
  if (runParallel) {
    jnk3<- future_pmap(list(combinedE$nxTo,combinedE$nyTo,combinedE$numTo),deDuplicateToDataLemma,.progress=FALSE)
  } else {
    jnk3<- pmap(list(combinedE$nxTo,combinedE$nyTo,combinedE$numTo),deDuplicateToDataLemma)
  }
  if (debugOutput){toc()}
  
  #now put the new nxTo,nyTo and numTo back into combinedE
  combinedE$nxTo<-map(jnk3,function(a){c(a[1],recursive=TRUE)})
  combinedE$nyTo<-map(jnk3,function(a){c(a[2],recursive=TRUE)})
  combinedE$numTo<-map(jnk3,function(a){c(a[3],recursive=TRUE)})
  
  #if either input had lat and lon data, add it back
  if (haveLatLon) {
    combinedE<-addLatLon(combinedE)
  }
  
  #done!
  return(combinedE)
}

#===============================================================================
# the following routines takes the connectivity data and uses it to calculate the
# dispersal from an arbitrary number of initial points with the function
# dispersePoints(). The points are stored in a data frame that has columns nx,
# ny and num. (nx,ny) are the location in the mercator grid of the particles,
# and num is the number of of individuals at those locations. A helper function
# addLatLon() adds the columns lat and lon to the data frame, to aid in
# plotting.

#this function takes a connectivity data structure and returns a dictionary 
#that maps from a (nx,ny) point to the row of the connectivity that has that
#(nx,ny) point as its (nxFrom,nyFrom) point
makeConnectionDict<-function(E){
  connectDict<-dict(1:nrow(E),map2(E$nxFrom,E$nyFrom,c)) #items, keys!
  return(connectDict)
}

#this function takes a data.frame orgDist with columns nx, ny and num that
#represent the distribution of organisms; the connectivity structure E; and a
#connection dictionary connectDict (made with makeConnectDict()). It returns a
#similar structure that contains the distribution of organisms one generation
#later, assuming a growth rate of R.
#
#The growth rate of R is the number of propagules release for a single 
#individual. It includes the effects of mortality while drifting, but DOES NOT
#include the effects of propagules leaving suitable habitat. 
propagateOneGeneration<-function(E,connectDict,orgDist,R) {
  
  #calculate the total number which when multiplied by the sum of numTo for each
  #row gives the number of particles released at each location in E
  #per individual at that location. Then multiply it by R
  RscaleFactor=R/(E$numLaunched[1]) #assumes numLaunched is the same for all rows!
  
  #now build the next generation. I am sure this could be done more efficiently...
  nextGen<-data.frame()
  for (nCrit in 1:nrow(orgDist)){
    #the as.integer is because c() by default makes num
    whichRow<-connectDict$get(as.integer(c(orgDist$nx[nCrit],orgDist$ny[nCrit])))
    newRecruits<-data.frame(nx=E$nxTo[[whichRow]],ny=E$nyTo[[whichRow]],num=E$numTo[[whichRow]]*orgDist$num[nCrit])
    nextGen<-rbind(nextGen,newRecruits)
  }
  #now, if there were multiple points, nextGen will have multiple nx,ny pairs. We need 
  #to collapse this by grouping by nx,ny pairs and adding their num's
  nextGen<-nextGen %>% group_by(nx,ny) %>% summarise(num=sum(num))
  
  #now normalize population to get total population correct
  nextGen$num<-nextGen$num*RscaleFactor
  
  return(nextGen)
}

#this function adds lat and lon to the orgDist data frame to make plotting
#easier. 
addLatLon2orgDist<-function(orgDist){
  returnVals<-nxny2lonlat(orgDist$nx,orgDist$ny)
  
  #do the "from", each a simple vector
  orgDist$lon<-returnVals$lonVec
  orgDist$lat<-returnVals$latVec
  
  return(orgDist)
}



#===============================================================================
#the code below here is for testing, and is wrapped in if-statements. Only 
#make the if's true to test the code above. 

#if True, test addLatLon
#this code also reads in the data used by subsequent test blocks
if (FALSE) {
  
  #get data to process
  regionName<-'theAmericas'
  depth<-1 #what depth the larvae are released from
  year<-'climatology' #the year you want to get, or "climatology" if you want all years
  verticalBehavior<-'fixed' #"fixed" if they stay at release depth, "starts" if they can advect vertically
  month<-5 #what month to get
  minPLD<-28; maxPLD<-minPLD #how many days they can drift. See text above for possible values
  
  #now get the data
  E<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
  
  Eplus<-addLatLon(E)
  if (FALSE) {
    #plot output
    plot(Eplus$lonFrom,Eplus$latFrom)
  }
}

#if True, test addLatLon AND the code to combine multiple connectivity files
#this code also reads in the data used by subsequent test blocks
if (FALSE) {
  
  #get data to process
  regionName<-'theAmericas'
  depth<-1 #what depth the larvae are released from
  year<-'climatology' #the year you want to get, or "climatology" if you want all years
  verticalBehavior<-'fixed' #"fixed" if they stay at release depth, "starts" if they can advect vertically
  month<-4 #what month to get
  minPLD<-28; maxPLD<-minPLD #how many days they can drift. See text above for possible values
  
  #now get the data for first month
  E1<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
  
  #and data for second month
  #IF MONTH IS UNCHANGED, numTo SHOULD DOUBLE WHEN MERGED... 
  month<-5 #KEEP AS SAME MONTH FOR SOME DEBUGGING, SEPERATE FOR OTHERS
  E2<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
  
  if (TRUE){
    #then subset data
    if (TRUE){
      #subset to a polygon
      lonLimits=c(-82.0,-80); latLimits=c(28.0,32.0)
      limitPoly=st_polygon(list(cbind(lonLimits[c(1,2,2,1,1)],latLimits[c(1,1,2,2,1)])))
      limitPoly<-st_sfc(limitPoly,crs=4326) #make into a surface, 4326 is WGS84
      E1<-addLatLon(E1)
      E2<-addLatLon(E2)
      E1<-subsetConnectivity_byPolygon(E1,limitPoly,trimTo=FALSE)
      E2<-subsetConnectivity_byPolygon(E2,limitPoly,trimTo=FALSE)
    } else {
      #subset to fixed number of lines, easy to debug...
      howMany<-1
      E1<-E1[1:howMany,]
      E2<-E2[1:howMany,]
    }
  }
  
  #combine data
  print('combining')
  E<-combineConnectivityData(E1,E2)
  print('done')
  
  Eplus<-addLatLon(E)
  if (FALSE) {
    #plot output
    plot(Eplus$lonFrom,Eplus$latFrom)
  }
}



#if true, test code to sub-set matrix by polygon, subsetConnectivity_byPolygon().
#to run requires test block that loads data to be run
if (FALSE) {
  #define limits to be included in connectivity
  #define as sf polygon, as described in https://gis.stackexchange.com/questions/403977/sf-create-polygon-from-minimum-x-and-y-coordinates
  lonLimits=c(-82.0,-80)
  latLimits=c(28.0,32.0)
  
  limitPoly=st_polygon(list(cbind(lonLimits[c(1,2,2,1,1)],latLimits[c(1,1,2,2,1)])))
  limitPoly<-st_sfc(limitPoly,crs=4326) #make into a surface, 4326 is WGS84
  
  #subset connectivity data
  EplusTrimmed<-subsetConnectivity_byPolygon(Eplus,limitPoly,trimTo=TRUE)
  
  #if true, try to plot limitPoly
  #this is a great resource
  if (TRUE) {
    library("rnaturalearth")
    library("rnaturalearthdata")
    
    world <- ne_countries(scale = "medium", returnclass = "sf")
    class(world)
     
    #plot(limitPoly)
    p<-ggplot(data = world) +
      geom_sf() +
      geom_sf(data=limitPoly,fill=NA) + #add limit polygon
      coord_sf(xlim= c(-88, -78), ylim = c(24.5, 33), expand = FALSE)
    #p<-p+geom_point(data = Eplus, aes(x = lonFrom, y = latFrom), size = 1, shape = 23, fill = "darkred")
    
    #now figure out how to add up all the points in lonTo and latTo, and count the
    #number of points to each lonTo and latTo point. Lets be sloppy, and iterate over rows
    allTo<-data.frame()
    for (n in 1:nrow(EplusTrimmed)){
      toFrame<-data.frame(lonTo=EplusTrimmed$lonTo[[n]],latTo=EplusTrimmed$latTo[[n]],
                          numTo=EplusTrimmed$numTo[[n]])
      allTo<-rbind(allTo,toFrame)
    }
    
    #now group and sum numTo
    whereAllWent<-allTo %>% group_by(lonTo,latTo) %>% summarise(numTo=sum(numTo))
    
    #now add plots 
    p<-p+geom_point(data = EplusTrimmed, aes(x = lonFrom, y = latFrom), size = 1, shape = 23, fill = "green")
    p<-p+geom_point(data = whereAllWent, aes(x = lonTo, y = latTo, fill=numTo,colour=numTo), size = 1, shape = 23)
    
    #now finish up plots
    print(p)
    
  }
}

#if true, test code to sub-set matrix by values on the model grid, subsetConnectivity_byGridValues().
#to run requires test block that loads data to be run
if (FALSE) {
  
  #define limits to be included in connectivity
  #define as sf polygon, as described in https://gis.stackexchange.com/questions/403977/sf-create-polygon-from-minimum-x-and-y-coordinates
  lonLimits=c(-87.0,-77.0)
  latLimits=c(28.0,30.0)
  
  lonLimits=c(-82.0,-80)
  latLimits=c(28.0,32.0)
  
  limitPoly=st_polygon(list(cbind(lonLimits[c(1,2,2,1,1)],latLimits[c(1,1,2,2,1)])))
  limitPoly<-st_sfc(limitPoly,crs=4326) #make into a surface, 4326 is WGS84
  
  #subset connectivity data to that within the a polygon
  EplusTrimmed<-subsetConnectivity_byPolygon(Eplus,limitPoly)  
  
  #subset connectivity data by grid value
  values<-gridDist; minVal<-2.2; maxVal<-10 #limit to distance from coast in grid cells
  values<-gridDepth; minVal<- 0; maxVal<-25 #limit to certain depth values
  EplusTrimmed<-subsetConnectivity_byGridValue(EplusTrimmed,values,minVal,maxVal,trimTo=TRUE)
  
  #get the depth data of model points to plot below. 
  theLons<-EplusTrimmed$lonFrom
  theLats<-EplusTrimmed$latFrom
  theDepths<-gridDepth[cbind(EplusTrimmed$nxFrom,EplusTrimmed$nyFrom)]
  if (TRUE){
    indx<-theDepths<6000
    theLons<-theLons[indx]
    theLats<-theLats[indx]
    theDepths<-theDepths[indx]
  }
  theData<-data.frame(theLons=theLons,theLats=theLats,theDepths=theDepths)
  
  #if true, plot results
  if (TRUE) {
    library("rnaturalearth")
    library("rnaturalearthdata")
    
    world <- ne_countries(scale = "medium", returnclass = "sf")
    class(world)
    
    #plot(limitPoly)
    p<-ggplot(data = world) +
      geom_sf() +
      geom_sf(data=limitPoly,fill=NA) + #add limit polygon
      coord_sf(xlim= c(-88, -78), ylim = c(24.5, 33), expand = FALSE)
    #p<-p+geom_point(data = Eplus, aes(x = lonFrom, y = latFrom), size = 1, shape = 23, fill = "darkred")

    #now figure out how to add up all the points in lonTo and latTo, and count the
    #number of points to each lonTo and latTo point. Lets be sloppy, and iterate over rows
    allTo<-data.frame()
    for (n in 1:nrow(EplusTrimmed)){
      toFrame<-data.frame(lonTo=EplusTrimmed$lonTo[[n]],latTo=EplusTrimmed$latTo[[n]],
                          numTo=EplusTrimmed$numTo[[n]])
      allTo<-rbind(allTo,toFrame)
    }
    #now group and sum numTo
    whereAllWent<-allTo %>% group_by(lonTo,latTo) %>% summarise(numTo=sum(numTo))
    
    #now add plots 
    p<-p+geom_point(data = EplusTrimmed, aes(x = lonFrom, y = latFrom), size = 1, shape = 23, fill = "green")
    p<-p+geom_point(data = whereAllWent, aes(x = lonTo, y = latTo, fill=numTo,colour=numTo), size = 1, shape = 23)
    #p<-p+geom_point(data = theData, aes(x = theLons, y = theLats, fill=theDepths,colour=theDepths), size = 1, shape = 23)
    #p<-p+geom_point(data = EplusTrimmed, aes(x = lonFrom, y = latFrom), size = 1, shape = 23, fill = "green")
    
    
    print(p)
    
  }
}

#this code debugs the progagateOneGeneration function
if (FALSE) {
  regionName<-'theAmericas'
  depth<-1
  year<-'climatology'
  verticalBehavior<-'fixed'
  month<-5
  minPLD<-18; maxPLD<-minPLD
  E<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
  Eplus<-addLatLon(E)
  #plot(Eplus$lonFrom,Eplus$latFrom)
  
  if (TRUE) {
    #complex test, multiple generations
    #now subset the data to a region and a depth
    #define as sf polygon, as described in https://gis.stackexchange.com/questions/403977/sf-create-polygon-from-minimum-x-and-y-coordinates
    lonLimits=c(-82.0,-79)
    latLimits=c(25.0,32.0)
    limitPoly=st_polygon(list(cbind(lonLimits[c(1,2,2,1,1)],latLimits[c(1,1,2,2,1)])))
    limitPoly<-st_sfc(limitPoly,crs=4326) #make into a surface, 4326 is WGS84
    
    #subset connectivity data to that within the a polygon
    EplusTrimmed<-subsetConnectivity_byPolygon(Eplus,limitPoly)  
    
    #subset connectivity data by grid value
    #it is important to include the trimTo=TRUE option, so all of the "To" points are also 
    #in the "From" points. 
    values<-gridDist; minVal<-2.2; maxVal<-10 #limit to distance from coast in grid cells
    values<-gridDepth; minVal<- 0; maxVal<-25 #limit to certain depth values
    EplusTrimmed<-subsetConnectivity_byGridValue(EplusTrimmed,values,minVal,maxVal,trimTo=TRUE)
    
    #now lets choose an initial point to seed the domain with. 
    orgDistInit<-data.frame(nx=2477,ny=1849,num=1)
    orgDistInit<-data.frame(nx=c(2477,2475),ny=c(1849,1891),num=c(0.5,0.5))
    
    #now lets see where it goes in two generations, using R of 10
    R<-1.0
    connectDict<-makeConnectionDict(EplusTrimmed)
    orgDistNext<-propagateOneGeneration(EplusTrimmed,connectDict,orgDistInit,R)
    orgDistNextNext<-propagateOneGeneration(EplusTrimmed,connectDict,orgDistNext,R)
    print(paste('The total population over three generations',
                sum(orgDistInit$num),sum(orgDistNext$num),sum(orgDistNextNext$num)))
    
    #now plot all the points in the habitat
    library("rnaturalearth")
    library("rnaturalearthdata")
    world <- ne_countries(scale = "medium", returnclass = "sf")
    class(world)
    
    #add latLon to orgDist data
    orgDistInit<-addLatLon2orgDist(orgDistInit)
    orgDistNext<-addLatLon2orgDist(orgDistNext)
    orgDistNextNext<-addLatLon2orgDist(orgDistNextNext)
    
    #plot(limitPoly)
    p<-ggplot(data = world) +
      geom_sf() +
      geom_sf(data=limitPoly,fill=NA) + #add limit polygon
      coord_sf(xlim= c(-88, -78), ylim = c(24.5, 33), expand = FALSE)
    p<-p+geom_point(data = EplusTrimmed, aes(x = lonFrom, y = latFrom), size = 1, shape = 23, fill = "darkred",
                    alpha=0.2)
    p<-p+geom_point(data = orgDistNextNext, aes(x = lon, y = lat), size = 1, shape = 23, fill = "darkseagreen1",
                    alpha=1.0)
    p<-p+geom_point(data = orgDistNext, aes(x = lon, y = lat), size = 1, shape = 23, fill = "firebrick1",
                    alpha=1.0)
    p<-p+geom_point(data = orgDistInit, aes(x = lon, y = lat), size = 1, shape = 23, fill = "green",
                    alpha=1.0)
    print(p)
  } else {
    #a super simple one generation sanity check
    #now lets choose an initial point to seed the domain with. 
    orgDistInit<-data.frame(nx=2477,ny=1849,num=1)
    orgDistInit<-data.frame(nx=c(2477,2475),ny=c(1849,1891),num=c(0.5,2.5))
    
    #now lets see where it goes in two generations, using R of 10
    R<-1.0
    connectDict<-makeConnectionDict(E)
    orgDistNext<-propagateOneGeneration(E,connectDict,orgDistInit,R)
    
    print(paste('total pop before and after for R',R,'is',sum(orgDistInit$num),sum(orgDistNext$num)))
  }
}