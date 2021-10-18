## ARI_IoU_Calculation

#set path were the simulation and results data will be
path <- "/Users/dan/Desktop/ARI_IoU_Indexes/test/"
#size of the square regions of interest to be analysed in nm. We recommend no larger than 2000x2000 nm.
ROIsize <- 2000
setwd(path)
#source file should be place in the same directory as "path"
source("ARI_IoU_source.R")
#name of the folders for each different simulation scenario
#pattern for this example is "simulation_"
#if you have named you folders differently please modify to match them.
pfolders<-list.files(path=path, pattern = "simulation_")

for (s in 1:length(pfolders)){
  simpath<- paste(path,paste(pfolders[s],"/",sep = ""), sep = "")
  simfolders<- mixedsort(list.files(simpath))
  for (s2 in 1:length(simfolders)){
    dpath<- paste(simpath,simfolders[s2], sep = "")
    cpath<- paste(simpath,paste(simfolders[s2],"/classes", sep = ""), sep = "")
    gtclass <- list.files(dpath, pattern="*.csv", full.names="TRUE")
    GTdata <- read.csv(gtclass)
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
    write.csv(cbind(class_names,ARI_values,IoU_values),paste(spath,paste(simfolders[s2], "ARI and IoU Results.csv", sep = "")))
    base::save.image(paste(spath,paste(simfolders[s2], "ARI and IoU Results.RData", sep = ""),sep = ""))
    #plots IoU polygons for each class
    polygon_plots(GTdata, IoU_output,spath,simfolders[s2],class_names)
  }
}