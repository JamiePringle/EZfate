# This code downloads connectivity data IF IT HAS NOT ALREADY BEEN DOWNLOADED.
# If you wish to download files which have been updated on the server, you 
# must delete them from the local directory first. 

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


getConnectivityData<-function(dataDir,regionName,depth,year,verticalBehavior,month,minPLD,maxPLD){
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
    download.file(dataURL,downloadToFile)
  }
  
  #load the data.frame and return it
  E<-readRDS(downloadToFile)
  return(E)
}


#====================================================================================
#the following is test code. Set conditional to false to use this as a module
if (FALSE) {
  #set up data for call
  dataDir<-'theTestData'
  regionName<-'theAmericas'
  depth<-40
  year<-'climatology'
  verticalBehavior<-'fixed'
  month<-2
  minPLD<-18; maxPLD<-18
  
  getConnectivityData(dataDir,regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
}