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
  # of points that went from (nxFrom,nyFrom) to (nxTo,nyTo)
  transDict=dict()
  
  #first we write the serial version of this code which iterates over rows in 
  #E
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
  
  #ok, transDict has the inverted matrix. Now lets invert it
  jnk<-Etrans$keys()
  jnk3<-to_list(for(p in jnk) as.numeric(p[[2]]))
  df<-data.frame(nxTo=I(jnk3))
  
  return(transDict)
}


#for debugging, lets run it with a given data set
if (TRUE) {
  regionName<-'theAmericas'
  depth<-1
  year<-2007
  verticalBehavior<-'starts'
  month<-5
  minPLD<-18; maxPLD<-minPLD
  E<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
  
  print('starting transpose')
  Etrans<-transposeConnectivity(E)
  print('done')
}