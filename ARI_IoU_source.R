## ARI and IoU source codes and libraries
library(mclustcomp)
library(dbscan)
library(gtools)
library(imager)
library(spatstat)
library(contoureR)
library(sp)
library(RColorBrewer)
library(grDevices)

ARI_calc <- function(GTdata,class){
  ARI_values <- c(rep_len(0, length.out = length(class)))
  for (i in 1:length(class)) {
    ARI<- mclustcomp::mclustcomp(GTdata$index,class[[i]]$result, type = "adjrand")
    ARI_values[i]<- ARI[1,2]
  }
  return(ARI_values)
} 

IoU_calc <- function(GTdata,class,mask_dim){
  IoU_values <- c(rep_len(0, length.out = length(class)))
for (ii in 1:length(class)){
    moleculeCoords <- cbind(GTdata$x, GTdata$y)
    moleculeClusterIndex<- GTdata$index
    #limits
    xmin <- min(moleculeCoords[,1])
    xmax <- max(moleculeCoords[,1])
    ymin <- min(moleculeCoords[,2])
    ymax <- max(moleculeCoords[,2])
    
    # preallocations and set up for GT loop
    NumInd <- unique(moleculeClusterIndex)
    NumClust <- length(NumInd[NumInd>0])
    ClustBoundInds <- vector("list", NumClust)
    boundptsGT <- vector("list", NumClust)
    GTpolys <- vector("list", NumClust)
    GTMasks <- vector("list", NumClust)
    GTmats <- vector("list", NumClust)
    
    # ground truth polygons
    for (u in 1:NumClust){
      cluster <- moleculeCoords[moleculeClusterIndex==u,1:2]
      ClustBoundInds[[u]] <- chull(cluster)
      
      if(length(ClustBoundInds[[u]])<3){
        GTmats[[u]] <- matrix(data = 0, nrow = mask_dim, ncol = mask_dim)
      }
      
      else{
        bpts <- cluster[ClustBoundInds[[u]],1:2]
        ACbptsInd <- orderPoints(bpts[,1], bpts[,2], xm = mean(bpts[,1]), ym = mean(bpts[,2]), clockwise = FALSE)
        boundptsGT[[u]] <-  list(x = bpts[ACbptsInd,1],y = bpts[ACbptsInd,2])
        GTpolys[[u]] <- owin(xrange=c(xmin,xmax), yrange=c(ymin,ymax), poly = boundptsGT[[u]])
        GTMasks[[u]] <- as.mask(GTpolys[[u]], eps=1, dimyx= mask_dim, xy=NULL)
        GTmats[[u]] <- GTMasks[[u]][["m"]]
      }
    }
    #plots GT mask (flattened)
    GT_full_mask <- Reduce('+', GTmats)
    flat_GT_mask <- GT_full_mask
    #plot(as.cimg(GT_full_mask))
    flat_GT_mask[which(GT_full_mask>0)]<-1
    #plot(as.cimg(flat_GT_mask))
    all_ResPolys<- vector("list", length(class))
    
    for (c in 1:length(class)){
      clusterRes<- class[[c]]$result
      counts <- table(clusterRes)
      clusters <- as.numeric(names(counts)[counts >= 3])[as.numeric(names(counts)[counts >= 3])>0]
      
      if (length(clusters) == 0){
        ResMats <- lapply(1, matrix, data = 0, nrow = mask_dim, ncol = mask_dim)
      }
      
      else{
        ResClust <- length(clusters)
        ResClustBoundsInds <- 0
        boundptsRes <- vector("list", ResClust)
        ResPolys <- vector("list", ResClust)
        ResMasks<- vector("list", ResClust)
        ResMats <- vector("list", ResClust)
        
        for (w in 1:length(clusters)){
          cluster <- moleculeCoords[which(clusterRes==clusters[w]),1:2]
          ResClustBoundsInds <- chull(cluster)
          bpts <- cluster[ResClustBoundsInds,1:2]
          ACbptsInd <- orderPoints(bpts[,1], bpts[,2], xm = mean(bpts[,1]), ym = mean(bpts[,2]), clockwise = FALSE)
          boundptsRes[[w]] <-  list(x = bpts[ACbptsInd,1],y = bpts[ACbptsInd,2])
          ResPolys[[w]] <- owin(xrange=c(xmin,xmax), yrange=c(ymin,ymax), poly = boundptsRes[[w]])
          ResMasks[[w]] <- as.mask(ResPolys[[w]], eps=1, dimyx=mask_dim, xy=NULL)
          ResMats[[w]] <- ResMasks[[w]][["m"]]
          #plot(as.cimg(ResMats[[w]]))
          all_ResPolys[[c]][[w]]<-ResPolys[[w]]
        }
      }
      
      Res_full_mask <- Reduce('+', ResMats)
      #plot(as.cimg(Res_full_mask))
      flat_Res_mask <- Res_full_mask
      #plot(as.cimg(Res_full_mask))
      flat_Res_mask[which(Res_full_mask>0)] <- 1
      #plot(as.cimg(flat_Res_mask))
      
      sum_mask <- flat_Res_mask+flat_GT_mask
      flat_union_mask <- sum_mask
      flat_union_mask[which(flat_union_mask>0)] <- 1
      overlap <- length(sum_mask[sum_mask==2])
      union <- length(flat_union_mask[flat_union_mask==1])
      Jac_IoU <- overlap/union
      IoU_values[c] <- Jac_IoU
    }
}
  IoU_output <- list(IoU_values = IoU_values, GTpolys = GTpolys, ResPolys = all_ResPolys)
  return(IoU_output)
}

polygon_plots <- function(GTdata, IoU_output, spath, FileName, class_names){
  pdf(paste(spath,paste(tools::file_path_sans_ext(FileName), "IoU Polygons.pdf", sep = " ")),width = 9, height = 10)
  magenta50 <- rgb(255, 0, 255, max = 255, alpha = 125, names = "magenta50")
  green50 <- rgb(0, 255, 0, max = 255, alpha = 125, names = "green50")
  dan_chrome <- c("#000000",rep(brewer.pal(12, "Paired"),50))
  clustcol <- dan_chrome[GTdata$index+1]
  for (p in 1:length(IoU_output$ResPolys)){
  par(mfrow=c(2,2))
  plot(cbind(GTdata$x,GTdata$y),col= clustcol, main = paste(tools::file_path_sans_ext(FileName), "Ground Truth Clusters"), xlab = "x (nm)", ylab = "y (nm)", pch = 19, cex= 0.5)
  #plots Ground Truth Polys
  plot(cbind(GTdata$x,GTdata$y),col= clustcol, main = "Ground Truth Polygons",xlab = "x (nm)", ylab = "y (nm)", pch = 19, cex= 0.5)
  lapply(IoU_output$GTpolys,function(x){plot.owin(x,add=TRUE,col=magenta50,border="magenta")})
  #plots Result Polys
  plot(cbind(GTdata$x,GTdata$y),col= clustcol, main = paste(tools::file_path_sans_ext(class_names[p]), "Class Result"), xlab = "x (nm)", ylab = "y (nm)", pch = 19, cex= 0.5)
  lapply(IoU_output$ResPolys[[p]],function(x){plot.owin(x,add=TRUE,col=green50,border="darkgreen")})
  #overlaid Polys for Result and GT
  plot(cbind(GTdata$x,GTdata$y),col= clustcol, main = "Overlay of Ground Truth and Result", xlab = "x (nm)", ylab = "y (nm)", pch = 19, cex= 0.5)
  lapply(IoU_output$ResPolys[[p]],function(x){plot.owin(x,add=TRUE,col=green50,border="darkgreen")})
  lapply(IoU_output$GTpolys,function(x){plot.owin(x,add=TRUE,col=magenta50,border="magenta")})
  }
  dev.off()
}
