# this module provides the function getFileFromOSN() which takes a
# file path as an arguement, and, if the data has not already been
# downloaded, downloads the corresponding file from the EZfate Open
# Storage Network S3 bucket, and stores it in the local directory
# OSNdataDir (or as specified immediately below. This directory will
# be created if it does not exist.
#
# If any of the files have been corrupted, or you just want to start
# over, you can delete OSNdataDir and all the files will be downloaded
# again when they are requested.

#set up the libraries we need
import s3fs
import os

#what directory to store data in
dataCacheDir='OSNdataDir/'

#information for EZfate OSN S3 storage
theBucket='ees230040-bucket01'
theToken='ULWSR4EC4O6GQB4225J2'
theEndpoint= "https://sdsc.osn.xsede.org"

def getFileFromOSN(fileName):
    '''getFileFromOSN(fileName) takes a single arguement, fileName, which
    is a string giving the path to a file in the EZfate data on the
    open storage network (OSN) S3 bucket. It then creates a matching
    path in the dataCacheDir cache directory (defined in the
    getEZfateFromOSN module, by default it is OSNdataDir), and
    downloads the file if it has not already been downloaded. The file
    will not be re-downloaded unless it is deleted from the cache directory.

    returns downloadFileTo, a string specifying the file to read to obtain 
    the data
    '''

    #check if path exists, if not, create path
    baseDir=os.path.dirname(dataCacheDir+fileName)
    if not os.path.exists(baseDir):
        print('directory',baseDir,'does not exists, and will be created')
        os.makedirs(baseDir)

    #check if file exists, if not, download file
    downloadFileTo=dataCacheDir+fileName
    if not os.path.exists(downloadFileTo):
        print('file',os.path.basename(fileName),'does not exist, downloading')
        fs=s3fs.S3FileSystem(anon=True,token=theToken,
                         client_kwargs={'endpoint_url': theEndpoint})
        #breakpoint()
        #print(fs.ls(theBucket+'/'+fileName))

        fs.download(theBucket+'/'+fileName,downloadFileTo)
        print('   done with %s, %d Mb downloaded'%(
            os.path.basename(fileName),os.path.getsize(dataCacheDir+fileName)/1024**2))
    else:
        print('file',os.path.basename(fileName),'exists, is in cache')

    return downloadFileTo

#now lets make sure that this works. Download one connectivity file and one grid file
if __name__ == "__main__":

    #first, lets get a grid file
    file2Download='EZfateData/EZfateFiles/ext-PSY4V3R1_mesh_zgr.nc'
    print('getting',file2Download)
    fileHere=getFileFromOSN(file2Download)
    print('   done, file located at',fileHere,'\n')


    file2Download=('EZfateData/communityConnectivityMatrices/theAmericas/1m/fixed/'+
                   'climatology_month03_minPLD30_maxPLD30.zip')
    print('getting',file2Download)
    fileHere=getFileFromOSN(file2Download)
    print('   done, file located at',fileHere,'\n')


    file2Download=('EZfateData/communityConnectivityMatrices/theAmericas/1m/fixed/'+
                   'climatology_month03_minPLD32_maxPLD32.zip')
    print('getting',file2Download)
    fileHere=getFileFromOSN(file2Download)
    print('   done, file located at',fileHere,'\n')
