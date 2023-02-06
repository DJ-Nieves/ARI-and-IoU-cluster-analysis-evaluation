### FOCAL clustering
library(gtools)
library(gplots)
library(dbscan)
library(raster)
library(imager)
library(spatstat)
library(pracma)
library(mclustcomp)

path <-"" #path to data files
setwd("") #output directory
dataFiles <- lapply(Sys.glob("*.csv"), read.csv)

for (f in 1:length(dataFiles)){
  library(imager)
  moleculeCoords <- cbind(dataFiles[[f]]$X, dataFiles[[f]]$X.1)
  moleculeClusterIndex <- dataFiles[[f]]$moleculeClusterIndex
  
  minL <- seq(1, 20, by = 1) 
  f_bins <- c((ceiling(max(moleculeCoords[,1])/10)+floor(abs(min(moleculeCoords[,1]))/10)),
              ceiling(max(moleculeCoords[,2])/10)+floor(abs(min(moleculeCoords[,1]))/10))
  
  xmin <- ceiling(min(moleculeCoords[,1]))
  xmax <- ceiling(max(moleculeCoords[,1]))
  ymin <- ceiling(min(moleculeCoords[,2]))
  ymax <- ceiling(max(moleculeCoords[,2]))
  
  histo <- hist2d(moleculeCoords[,1], moleculeCoords[,2],f_bins)
  counts <- histo[["counts"]]
  
  foc_mat <- matrix(data = 0, nrow = f_bins[1]+1, ncol = f_bins[2]+1)
  
  for (i in 2:length(counts[,1])-1){
    for (j in 2:length(counts[1,])-1){
      foc_mat[i+1,j+1] <- sum(counts[i,j],counts[i,j-1],counts[i,j+1],
                              counts[i-1,j],counts[i-1,j-1], counts[i-1,j+1],
                              counts[i+1,j], counts[i+1,j-1], counts[i+1,j+1])
    }
  }
  foc_im <- as.cimg(foc_mat)
  plot(foc_im)
  rand_scan2 <- matrix(0, nrow = length(minL))
  test_scan2 <- vector("list", length(minL))
  fmasks_filt2 <- vector("list", length(minL))
  for (r in 1:length(minL)){
  bordpts <- which(foc_mat<minL[r] & foc_mat>0, arr.ind = TRUE)
  corepts <- which(foc_mat>=minL[r], arr.ind = TRUE)
  
  if (length(corepts[,1])==0){
    f_mask <- matrix(data = 0, nrow = f_bins[1]+1, ncol = f_bins[2]+1)
  }
  
  if (length(corepts[,1])>0){
  if (length(bordpts[,1])==0){
    f_mask <- matrix(data = 0, nrow = f_bins[1]+1, ncol = f_bins[2]+1)
    if (length(corepts[,1])>0){
    for (c in 1:length(corepts[,1])){
      f_mask[corepts[c,1],corepts[c,2]] <-1
    }
    }
  }
  if (length(bordpts[,1])>0){
    bord_ind <- matrix(data=0, nrow = length(bordpts[,1]))
      for (b in 1:length(bordpts[,1])){
      brows<-c(bordpts[b,1]-1,bordpts[b,1]+1)
      bcols<-c(bordpts[b,2]-1,bordpts[b,2]+1)
      bord_reg<- foc_mat[brows[1]:brows[2],bcols[1]:bcols[2]]
      bord_ind[b]<- length(which(bord_reg>minL[r]))
      }
    f_mask <- matrix(data = 0, nrow = f_bins[1]+1, ncol = f_bins[2]+1)
    t_bord <- bordpts[bord_ind>0,1:2] #need if statment here to catch no borders
    if (length(t_bord)>0){ #new if here
      for (c in 1:length(corepts[,1])){
        f_mask[corepts[c,1],corepts[c,2]] <-1
      }
      for (tb in 1:length(t_bord[,1])){
        f_mask[t_bord[tb,1],t_bord[tb,2]] <-1
      }}
    if (length(t_bord)==0){ #this is new loop
      for (c in 1:length(corepts[,1])){
        f_mask[corepts[c,1],corepts[c,2]] <-1
      }
      }}
    }
  #plot(as.cimg(f_mask))
  range <- range(as.cimg(f_mask))
  if (range[1]+range[2]==0){
    F_clust <- rep(0, length(moleculeCoords[,1]))
    fmasks_filt2[[r]]<-NA
  }
  if (range[1]+range[2]==1){
    plot(as.cimg(flipdim(rot90(rot90(f_mask)))))
    fmasks <- split_connected(as.cimg(f_mask))
    filt <- matrix(0, length(fmasks))
    for (l in 1:length(fmasks)){
      filt[[l]] <- length(which(fmasks[[l]]==TRUE))
      fmasks_filt<- which(filt>9)
      fmasks_filt2[[r]]<-fmasks_filt
    }
    if(length(fmasks_filt==0)){
      F_clust <- rep(0, length(moleculeCoords[,1]))
    }
    if(length(fmasks_filt>0)){
      F_clust <- vector(mode = "integer",length = length(moleculeCoords[,1]))
      f_pts_in <- vector("list", length(fmasks_filt))
      for (k in 1:length(fmasks_filt)){
        mask <- owin(xrange = c(0,xmax+abs(xmin)), yrange = c(0,ymax+abs(ymin)),mask = rot90(flipdim(as.matrix(fmasks[[fmasks_filt[[k]]]]), dim = 2)))
        f_pts_in[[k]] <- inside.owin(moleculeCoords[,1]+abs(xmin), moleculeCoords[,2]+abs(ymin), mask)
        F_clust[which(f_pts_in[[k]]==TRUE)]<-k
      }}
  }
  test_scan2[[r]] <- F_clust
  rand_scan<- mclustcomp::mclustcomp(moleculeClusterIndex,test_scan2[[r]], type = "adjrand")
  rand_scan2[r] <- rand_scan[1,2]
  }
  fname <- paste(f,".RData", sep = "")
  detach("package:imager", unload=TRUE)
  save.image(file = paste(path,fname,sep = "")) #need to sort our saving of environment
  rm(list= ls()[!(ls() %in% c('dataFiles','path'))])
}
  