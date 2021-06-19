## ARI and IoU calculation for SMLM clustering results

# set dpath to the directory that your ground truth data (in the 3 column form: x,y,index) is located
# **NOTE** only one ground truth data file can be placed in the directory
path <- "/Users/dan/Desktop/ARI_IoU_Indexes/"
setwd(path)
#source file should be place in the same directory as the ground truth data file
source("ARI_IoU_source.R") 
#input name of folder containing the ground truth data file, and the name of the file with extension (.csv)
FolderName <- "test/"
dpath <- paste(path,FolderName, sep = "")
FileName <- "data.csv"
#size of the square regions of interest to be analysed in nm. We recommend no larger than 2000x2000 nm.
ROIsize <- 2000

#reads in ground truth data
GTdata <- read.csv(paste(dpath,FileName,sep = ""))

# all your cluster class results (indexes from clustering analysis for each point, single column) 
# should be place in a folder called "classes" within the ground truth data folder
# **NOTE** single .csv file for each cluster result.
cpath <- paste(dpath,"/classes", sep = "")
cfiles <- mixedsort(list.files(cpath, pattern="*.csv", full.names="TRUE"))
class <- lapply(Sys.glob(cfiles), read.csv)

#calculates Adjusted Rand Index for each class file
ARI_values <-ARI_calc(GTdata,class)

#calculates Intersection over Union for each class file
mask_dim <- ROIsize+400
IoU_output <- IoU_calc(GTdata,class, mask_dim)
IoU_values <- IoU_output$IoU_values
class_names <- basename(cfiles[1:length(cfiles)])

#creates new directory for ARI and IoU results to be saved
spath <- paste(cpath,paste("/ARI and IoU Result ",paste(format(Sys.time(), "%F %H_%M"),"/",sep = ""), sep = ""),sep = "")
dir.create(spath,showWarnings = TRUE, recursive = FALSE)
#writes and saves ARI and IoU results to .csv file
write.csv(cbind(class_names,ARI_values,IoU_values),paste(spath,paste(tools::file_path_sans_ext(FileName), "ARI and IoU Results.csv", sep = " "),sep = ""))
base::save.image(paste(spath,paste(tools::file_path_sans_ext(FileName), "ARI and IoU Results.RData", sep = " "),sep = ""))
#plots IoU polygons for each class
polygon_plots(GTdata, IoU_output, spath, FileName,class_names)