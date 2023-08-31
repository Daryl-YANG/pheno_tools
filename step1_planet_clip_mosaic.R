###########################################################################################
#
#  this script extracts spectral reflectance (and quality control layer if preferred) from
#  original modis brdf-corrected reflectance hdf files downloaded from DACC, convert hdf 
#  to tif and merge spectral bands. the output is a single raster file for each hdf file
#  that contains all desired spectral bands and quality control layers
#
#    --- Last updated:  2021.06.22 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("ggplot2","raster", 'rgdal', 'stringr')  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#************************************ user parameters ************************************#
# define an output directory to store outputs
outDIR <- file.path("G:\\My Drive\\github\\pheno_tools\\planetscope\\test")
# create output directory if not exist
if (! file.exists(outDIR)) dir.create(outDIR,recursive=TRUE)
# create an temporary to store files temporarily generated during the course of processing
tempDIR <- file.path(paste0(outDIR, "/", 'temporary'))
if (! file.exists(tempDIR)) dir.create(tempDIR,recursive=TRUE)

# define output project system
outTRS <- '+proj=utm +zone=03 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

# define string pattern to remove from the default planetscope file format
str_rm_pattern <- '_psscene_analytic_sr_udm2.zip'
#*****************************************************************************************#

#*************************************** load data ***************************************#
# define the directory to hdf files
dataDIR <- file.path("Z:\\dyang\\projects\\ngee_arctic\\phase4_data\\planetscope\\toolik\\2020")
# search all hdf files stored in the directory
zipLIST <- list.files(dataDIR, pattern = '*.zip', recursive = FALSE, all.files = TRUE, 
                      full.names = TRUE)
# check the number of hdf files
print(paste0('number of zip to process: ', length(zipLIST)))
# print out the first a couple of hdf files to check
str(zipLIST)

# load in shp file that defines the region to clip
# best only have one feature in the shp file
shpDIR <- file.path("G:\\My Drive\\projects\\ngee_arctic\\phase4\\evaluation_sites\\tier1\\box_tolik.shp")
shpVCT <- readOGR(shpDIR)
# project the vector file to target projection system to match with raster files
shpVCT <- spTransform(shpVCT, outTRS)
#*****************************************************************************************#

#*************************************** unzip file **************************************#
for (zip in zipLIST)
{
  # create a directory to store unzipped files
  foldername <- str_remove(basename(zip), str_rm_pattern)
  
  extrDIR <- file.path(outDIR, foldername, '/')
  if (! file.exists(extrDIR)) dir.create(extrDIR,recursive=TRUE)
  # extract site and date of data collection for output folder name
  unzip(zip, exdir = tempDIR, overwrite = TRUE)
  # check if there are multiple images collected on the data
  files <- list.files(file.path(tempDIR, 'PSScene', '/'), pattern = '.tif', 
                      recursive = TRUE, all.files = TRUE, 
                      full.names = TRUE)
  #numSR <- length(sr_reps)
  # group images if multiple images are collected on a same day
  dataID <- unique(substr(basename(files), 1, 15))
  for (id in dataID)
  {
    sub_outDIR <- file.path(extrDIR,id)
    dir.create(sub_outDIR)
    sub_fileLIST <- list.files(file.path(tempDIR, 'PSScene', '/'), pattern = id,
                               recursive = TRUE, all.files = TRUE, 
                               full.names = TRUE)
    file.copy(sub_fileLIST,sub_outDIR)
  }
  ### delete all files in the temporary folder, so that it can be used for next 
  ### round hdf processing
  unlink(tempDIR, recursive = T, force = T)
}
#*****************************************************************************************#

#************************************** mosaic tiles *************************************#
unzippedLIST <- dir(outDIR, full.names = TRUE)
for(folder in unzippedLIST)
{
  # search all surface reflectance files in the folder
  sr_files <- list.files(folder, pattern = '_SR_', recursive = TRUE, 
                         all.files = TRUE, full.names = TRUE)
  # if surface reflectance files are found, then mosaic them
  if (length(sr_files) == 1)
  {
    mosaicRST <- brick(sr_files)
    mosaicRST <- projectRaster(mosaicRST, crs=outTRS)
    mosaicRST <- crop(mosaicRST, shpVCT, snap = 'out')
    # search quality layer and apply mask
    qaDIR <- gsub('AnalyticMS_SR_harmonized', 'udm2', sr_files)
    mosaic.qaRST <- brick(qaDIR)
    mosaic.qaRST <- projectRaster(mosaic.qaRST, crs=outTRS)
    mosaic.qaRST <- crop(mosaic.qaRST, shpVCT, snap = 'out')
    # create mask
    shadow <- mosaic.qaRST[[3]]
    hazelight <-  mosaic.qaRST[[4]]
    hazeheavy <- mosaic.qaRST[[5]]
    cloud <- mosaic.qaRST[[6]]
    mask.INT <- shadow + hazelight + hazeheavy + cloud
    mask.INT[mask.INT > 0] <- 1
    # apply the mask
    mosaicRST.masked <- mask(mosaicRST, mask.INT, maskvalue = 1)
    
    foldername <- basename(folder)
    outNAME <- paste0(folder, '/',  foldername, '_mosaicqc.tif')
    writeRaster(mosaicRST.masked, outNAME, format="GTiff", overwrite=TRUE)
  }
  
  if (length(sr_files) > 1)
  {
    sr_RSTs <- list()
    for(i in 1:length(sr_files)) { sr_RSTs[i] <- brick(sr_files[i]) }
    
    sr_RSTs$fun <- mean
    mosaicRST <- do.call(mosaic,sr_RSTs)
    # convert mosaic file to target projection system
    mosaicRST <- projectRaster(mosaicRST, crs=outTRS)
    mosaicRST <- crop(mosaicRST, shpVCT, snap = 'out')
    
    # search quality layer and apply mask
    qaDIRs <- gsub('AnalyticMS_SR_harmonized', 'udm2', sr_files)
    qa_RSTs <- list()
    for(i in 1:length(qaDIRs)) { qa_RSTs[i] <- brick(qaDIRs[i]) }
    qa_RSTs$fun <- mean
    mosaic.qaRST <- do.call(mosaic, qa_RSTs)
    mosaic.qaRST <- projectRaster(mosaic.qaRST, crs=outTRS)
    mosaic.qaRST <- crop(mosaic.qaRST, shpVCT, snap = 'out')
    #mosaic.qaRST[mosaic.qaRST == 0] <- NA
    
    # create mask
    shadow <- mosaic.qaRST[[3]]
    hazelight <-  mosaic.qaRST[[4]]
    hazeheavy <- mosaic.qaRST[[5]]
    cloud <- mosaic.qaRST[[6]]
    mask.INT <- shadow + hazelight + hazeheavy + cloud
    mask.INT[mask.INT > 0] <- 1
  
    # apply the mask
    mosaicRST.masked <- mask(mosaicRST, mask.INT, maskvalue = 1)
    
    # save mosaic raster
    foldername <- basename(folder)
    outNAME <- paste0(folder, '/',  foldername, '_mosaicqc.tif')
    writeRaster(mosaicRST.masked, outNAME, format="GTiff", overwrite=TRUE)
  }
}
#*****************************************************************************************#







