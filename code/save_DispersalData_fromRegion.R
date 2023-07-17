source('connectivityUtilities.R')
library(sf)
library(ggplot2)
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf") #get coastline data
class(world)

#make a data frame that contains where all points went as a function of PLD for a 
#given polygon. The resulting data frame will contain this data for many PLD, and 
#the data frame will have each row labled with the PLD, so we can animate it below. This 
#will be called whereAllWentAllPLD
whereAllWentAllPLD=data.frame()

for (minPLD in seq(2,18,2)) {
  tic(paste('working on PLD',minPLD))
  print(paste('starting PLD',minPLD))
  Eall<-data.frame()
  tic('getting all the data')
  for (regionName in c('theAmericas')) {
    print(paste('working on',regionName))
    
    #set up geographical subregion to get data from
    #define polyName, which specifies the root of the filename
    #the data will be stored in 
    if (FALSE) {
      polyName<-'EastFlorida'
      lonLimits=c(-82.0,-80) #longitude of corners
      latLimits=c(28.0,30.0) #latitude of corners
      limitPoly=st_polygon(list(cbind(lonLimits[c(1,2,2,1,1)],latLimits[c(1,1,2,2,1)])))
    } else if (FALSE) {
      polyName<-'QuintanaRoo' #from 17.491481, -88.148592
      lonLimits=c(-88.14,-88.14+2.0) #longitude of corners
      latLimits=c(17.49,17.49+2.0) #latitude of corners
      limitPoly=st_polygon(list(cbind(lonLimits[c(1,2,2,1,1)],latLimits[c(1,1,2,2,1)])))
    } else if (TRUE) {
      polyName<-'CODE_region' #Coastal Ocean Dynamics Region, 
      #coastal points 38.27, -122.99 and 40.29, -124.35
      lonLimits=c(-123.0-2.0,123.0) #longitude of corners
      latLimits=c(38.272,38.272+2.0) #latitude of corners
      limitPoly=st_polygon(list(cbind(c(-122.99,-124.25,-126.35,-124.99,-122.99),
                                      c(  38.27,  40.29,  40.29,  38.27,  38.27))))
    }
    
    #define depth and year/climatology and vertical behavior. This information
    #will be added to polyName
    if (TRUE) {
      depth<-10; depthName=paste('_',as.character(depth),'m',sep='')
      year<-'climatology'
      verticalBehavior<-'fixed' #'starts'
    }
    
    #add vertical behavior information to polyName
    polyName<-paste(polyName,year,depthName,verticalBehavior,sep='_')
    
    #get data
    tic('getting data')
    if (TRUE) {
      timeName='AprMayJune_'
      maxPLD<-minPLD
      month<-4; E1<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
      month<-5; E2<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
      month<-6; E3<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
    } 
    toc()
    
    #trim to land
    #tic(paste('have trimmed',regionName))
    #trimTo=TRUE
    #E1<-subsetConnectivity_byGridValue(E1,gridDist,0.0,2.1,trimTo=trimTo)
    #E2<-subsetConnectivity_byGridValue(E2,gridDist,0.0,2.1,trimTo=trimTo)
    #E3<-subsetConnectivity_byGridValue(E3,gridDist,0.0,2.1,trimTo=trimTo)
    #toc() 
    
    #now add lat and lon to all the files, and then subset to small region
    #define limits of a box that defines the release locations
    tic('add lat lons')
    limitPoly<-st_sfc(limitPoly,crs=4326) #make into a surface, 4326 is WGS84
    E1<-addLatLon(E1)
    E2<-addLatLon(E2)
    E3<-addLatLon(E3)
    toc()
    
    #now make file name for output
    dataFileName<-paste('connectivityData/',timeName,polyName,sep='')
    
    tic('subset by polygon')
    E1<-subsetConnectivity_byPolygon(E1,limitPoly,trimTo=FALSE)
    E2<-subsetConnectivity_byPolygon(E2,limitPoly,trimTo=FALSE)
    E3<-subsetConnectivity_byPolygon(E3,limitPoly,trimTo=FALSE)
    toc()
    
    
    tic(paste('have combined',regionName))
    E<-combineConnectivityData(E1,E2)
    E<-combineConnectivityData(E,E3)
    toc()
    
    if (nrow(Eall)==0) {
      Eall<-E
    }
    else {
      tic(paste('have combined',regionName,'with other regions'))
      Eall<-combineConnectivityData(Eall,E)
      toc()
    }
  }
  toc()
  
  tic('making allTo and whereAllWent')
  #now lets combine each row in Eall and make data.frame with each location particles get
  #to and the total number of particles that get there.
  allTo<-data.frame(lonTo=c(Eall$lonTo,recursive=TRUE),
                    latTo=c(Eall$latTo,recursive=TRUE),
                    numTo=c(Eall$numTo,recursive=TRUE))
  
  #and them add together the total number of points that get to the same (lonTo,latTo) points
  whereAllWent<-allTo %>% group_by(lonTo,latTo) %>% summarise(numTo=sum(numTo))
  
  #normalize by number that came to this point
  whereAllWent$fracTo<-whereAllWent$numTo/sum(Eall$numLaunched)
  
  #add a PLD column
  whereAllWent$PLD<-minPLD
  toc()
  
  tic(paste('adding PLD',minPLD,' to whereAllWentAllPLD'))
  #and add to whereAllWentAllPLD
  whereAllWentAllPLD<-rbind(whereAllWentAllPLD,whereAllWent[c('lonTo','latTo',
                                                              'numTo','fracTo','PLD')])
  toc()
  
  #let people know how long things are taking per PLD
  toc()
}

#{r save data to RDS}
saveRDS(whereAllWentAllPLD,file=paste(dataFileName,'.RDS',sep=''))

#now make a unique data frame with all starting locations of last time
Estart=E[c('lonFrom','latFrom')]
saveRDS(Estart,file=paste(dataFileName,'_startPoints.RDS',sep=''))
