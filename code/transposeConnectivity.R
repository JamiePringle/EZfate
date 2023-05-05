#===============================================================================
# This code defines a function transposeConnectivity() that takes a connectivity
# matrix and creates its transpose. In this transposed matrix, nxTo and nyTo are
# the model grid cells where the Lagrangian paths end, and in the data.frame()
# are a list of integers, as are nxFrom and nyFrom in the original matrix E. The
# columns nxFrom and nyFrom in Etrans are now lists of all the locations the
# pathways came from, similar to nxTo and nyTo in the original matrix E.
#
# This code also provides a helper function, addLatLon2transpose() that does
# what addLatLon does to a forward in time matrix, but to the transposed matrix.

# THIS CODE ASSUMES connectivityUtilities has already been sourced, so that
# all the libraries it loads are present. 

library(comprehenr)

transposeConnectivity<-function(E) {
  
  # the first step is to make a dictionary that has as keys all (nxTo,nyTo) pairs
  # in E, and as values a dictionary. The keys to that second dictionary are a
  # (nxFrom,nyFrom) pair that went to (nxTo,nyTo), and the value are the number 
  # of points that went from (nxFrom,nyFrom) to (nxTo,nyTo) (like numTo, but numFrom)
  transDict=dict()
  
  #first we write the serial version of this code which iterates over rows in 
  #E. The point of this code is to be clear and obvious. I can work on speed
  #once it is correct. Knuth et al. 
  for (n in 1:nrow(E)){
    theRow<-E[n,]
    
    #this is the key to the inner dictionary
    innerKey=c(theRow$nxFrom,theRow$nyFrom)
    
    #now loop over (nxTo,nyTo) pairs. Create a new entry into transDict if that
    #key does not exist their, or update existing entry if it does. Assumes, as 
    #should be true, that nxTo,nyTo and numTo have same length
    for (nn in 1:length(theRow$numTo)){
      outerKey=c(theRow$nxTo[[1]][nn],theRow$nyTo[[1]][nn])
      thisDict=transDict$get(outerKey,dict()) #get dict
      thisDict$set(innerKey,thisDict$get(innerKey,0)+theRow$numTo[[1]][nn]) #modify
      transDict$set(outerKey,thisDict) #set it again
    }
  }
  print('done making dictionary, now make data.frame')
  
  #I AM ASSUMING THAT THE ORDER OF ITEMS RETURNED FROM DICT IS STABLE
  #OVER TIME IF I DO NOT INSERT ANYTHING INTO IT BETWEEN ACCESS...
  
  #ok, transDict has the inverted matrix. Now lets invert it
  transKeys<-transDict$keys()
  #transValues<-transDict$values()
  nxTo<-to_list(for(p in transKeys) as.numeric(p[[1]]))
  nyTo<-to_list(for(p in transKeys) as.numeric(p[[2]]))
  Etrans<-data.frame(nxTo=I(nxTo),nyTo=I(nyTo))
  
  #now make columns of empty lists for remaining fields.
  #Etrans$nxFrom<-rep(list(list()),nrow(Etrans))
  #Etrans$nyFrom<-rep(list(list()),nrow(Etrans))  
  #Etrans$numFrom<-rep(list(list()),nrow(Etrans))
  Etrans$nxFrom<-rep(c(-1),nrow(Etrans))
  Etrans$nyFrom<-rep(c(-1),nrow(Etrans))  
  Etrans$numFrom<-rep(c(-1),nrow(Etrans))
  
  class(Etrans$nxFrom)<-"vector"
  class(Etrans$nyFrom)<-"vector"
  class(Etrans$numFrom)<-"vector"
  
  #now loop over transDict, and fill in lists. I know this is slow...
  #but I can figure out its correctness easilly
  for (n in 1:length(transKeys)) {
    theKey<-transKeys[[n]]
    innerDict<-transDict$get(theKey)
    innerKeys<-innerDict$keys()
    innerValues<-innerDict$values()
    nxFromList<-to_list(for(p in innerKeys) as.numeric(p[[1]]))
    nyFromList<-to_list(for(p in innerKeys) as.numeric(p[[2]]))
    numFromList<-to_list(for(p in innerValues) as.numeric(p[[1]]))
    
    #no, I have no idea exactly why this works to put a list
    #into a signle cell of a data.frame()
    #Etrans[n,'nxFrom'][[1]]<-I(list(c(nxFromList,recursive=TRUE)))
    #Etrans[n,'nyFrom'][[1]]<-I(list(c(nyFromList,recursive=TRUE)))
    #Etrans[n,'numFrom'][[1]]<-I(list(c(numFromList,recursive=TRUE)))
    
    Etrans$nxFrom[n]<-list(c(nxFromList,recursive=TRUE))
    Etrans$nyFrom[n]<-list(c(nyFromList,recursive=TRUE))
    Etrans$numFrom[n]<-list(c(numFromList,recursive=TRUE))
  }
  
  #the nxTo and nyTo should be just numbers, not lists, to match E not transposed
  class(Etrans$nxTo)<-'numeric'
  class(Etrans$nyTo)<-'numeric'
  
  return(Etrans)
}

#this function takes a TRANSPOSED dispersal structure and returns its with all nx* and ny* variables 
#suplemented with lon* and lat* values
addLatLon2transpose<-function(E){
  #nxny2lonlat is defined in connectivityUtilities
  returnVals<-nxny2lonlat(E$nxTo,E$nyTo)
  
  #do the "from", each a simple vector
  E$lonTo<-returnVals$lonVec
  E$latTo<-returnVals$latVec
  
  # #apply to each row for lonTo and latTo
  # print('start apply 1')
  # E$lonTo<-apply(E,1,nxny2lon)
  # print('end apply 2')
  # E$latTo<-apply(E,1,nxny2lat)
  # print('done with applies')
  
  #print('hello')
  
  # this is not memory efficient, but is it fast?
  nxFromList<-apply(E,1,function(a) {return(a$nxFrom)})
  nyFromList<-apply(E,1,function(a) {return(a$nyFrom)})
  E$lonFrom<-map2(nxFromList,nyFromList,nxny2lon)
  E$latFrom<-map2(nxFromList,nyFromList,nxny2lat)
  
  return(E)
}

#for debugging, lets run it with a given data set
if (FALSE) {
  if (FALSE) {
    #load from raw data file
    regionName<-'theAmericas'
    depth<-1
    year<-'climatology'
    verticalBehavior<-'starts'
    month<-5
    minPLD<-18; maxPLD<-minPLD
    E<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
  } else {
    #use precomputed E
    E<-readRDS('mayJuneEastCoast_connectivity.RDS')
  }
  print('starting transpose')
  Etrans<-transposeConnectivity(E)
  print('done, now add lat and lon')
  Etrans<-addLatLon2transpose(Etrans)
}
#===============================================================================