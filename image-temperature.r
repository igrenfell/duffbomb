library(EBImage)
library(lattice)
library(viridisLite)
library(RColorBrewer)
library(raster)
library(rasterVis)
setwd("I:\\workspace\\smouldering")

setwd("I:\\workspace\\smouldering\\biochar2-06252024")
setwd("I:\\workspace\\smouldering\\biochar3")
setwd("I:\\workspace\\smouldering\\biochar5")
setwd("I:\\workspace\\smouldering\\biochar4")
setwd("I:\\workspace\\smouldering\\october-images")
setwd("I:\\workspace\\smouldering\\nitrogen")


allpngs <- Sys.glob("*.png")
alltifs <- Sys.glob("*.tiff")
# 
# rootstr <- "charcoal1sm"
# 
# flist <- allpngs[grep(rootstr, allpngs)]
flist <- alltifs

nfiles <- length(flist)


#y <-  55.621x2 - 590.08x + 2193.2
#y <- 56.552x2 - 584.38x + 2177.7



rotate <- function(x) t(apply(x, 2, rev))
coul <- viridis(100)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
tempvec <- rep(NA, nfiles)
medvec <- rep(NA, nfiles)
maxvec <- rep(NA, nfiles)

for(i in 1:nfiles)
{
  
  setwd("I:\\workspace\\smouldering\\nitrogen")
  
  
  f <- flist[i]
  
  fsplit <- strsplit(f, "-")[[1]]
  wspeed <- fsplit[6]
  wspeed <- gsub(".tiff", "", wspeed)
  
  img <- readImage(f)
  rband <- img[,,1]
  gband <- img[,,2]
  bband <- img[,,3]
  
  rvec <- as.numeric(rband)
  bvec <- as.numeric(bband)
  gvec <- as.numeric(gband)
  
  svec <- rvec + gvec
  ratband <- rband / gband
  
  rdat <- imageData(rband)
  bdat <- imageData(bband)
  gdat <- imageData(gband)
  
  rraster <- raster(rdat)
  braster <- raster(rdat)
  graster <- raster(gdat)
  
  
  ratdat <- imageData(ratband)
  ratraster <- raster(rraster)
  ratvec <- as.numeric(ratdat)
  ratvec[is.infinite(ratvec)] <- .001
  ratvec[is.na(ratvec)] <- .001
  
  #ratvec <- na.exclude(ratvec)
  
  #transvec <- 55.621*ratvec**2 - 590.08*ratvec + 2193.2
  transvec <- 56.552*ratvec**2 -  584.38*ratvec +  2177.7
  transvec[svec < 0.0001] <- 0
  maxval <- 56.552 -  584.38 +  2177.7
  nnval <- quantile(transvec, 0.99)
  
  transvec[transvec > maxval] <- NA
  transsub <- transvec[transvec > 250]
  #transvec[bvec < .01 ] <- 0
  transvec[svec < .075 ] <- 0
  # transvec[transvec <1700] <- 0
  transvec[transvec <250] <- 0
  
  
  transmat <- matrix(transvec, nrow = dim(ratdat)[1], byrow = F)
  #transmat <- t(transmat)
  #transmat <- rotate(transmat)
  #transmat <- rotate(transmat)
  #transmat <- rotate(transmat)
  
  #transmat <- apply(transmat, 1, "rev")
  #transmat <- apply(transmat, 2, "rev")
  
  transimg <- Image(transmat)
  # plot(transimg)
  transraster <- raster(transmat)
  
  transsub <- transvec[transvec > 0]
  tempvec[i] <- mean(na.exclude(transsub))
  medvec[i] <- median(na.exclude(transsub))
  maxvec[i] <- max(na.exclude(transsub))
  subraster <- transraster
  subraster[subraster < 500] <- NA
  subraster <- t(subraster)
  transraster <- transraster[transraster]
  q95 <- quantile(na.exclude(transsub), 0.95)
  setwd("I:\\workspace\\smouldering\\nitrogen\\output")
  
  fout <- gsub(".tiff", "-contour.png", f)
  
  maxval <- 1500
  minval <- 1100
  my.at <- pretty(seq(1100, max(na.exclude(maxval)), length  = 6))
  
  png(fout, width = 2704 , height = 2028 , res = 300)
  print(levelplot(subraster, layers = 1,  zlim = c(1100, maxval), cex = 3, 
                  contour=TRUE, at = my.at, main = paste("95th Percentile  = ",
                  round(q95), ", wind = ", wspeed, sep = ""), margin = FALSE))
  dev.off()

  fout <- gsub(".tiff", "-histogram.png", f)
  png(fout)
  hist(na.exclude(transsub), breaks = 100, xlab = f, main = "", xlim = c(minval, maxval))
  dev.off()

  # is95 <- subraster
  # q95 <- quantile(na.exclude(transsub), 0.95)
  # is95[subraster >= q95] <- 1
  # is95[subraster < q95] <- 0
  # fout <- gsub(".tiff", "-95perc.png", f)
  # png(fout)
  # plot(is95, col = coul, main = paste("95th Percentile  = ", round(q95), ", wind = ", wspeed, sep = ""))
  # dev.off()
  # 
  # 
  # is1600 <- subraster
  # is1600[subraster >= 1600] <- 1
  # is1600[subraster < 1600] <- 0
  # fout <- gsub(".tiff", "-95perc.png", f)
  # png(fout)
  # plot(is1600, col = coul, main = paste("Temp > 1600 , wind = ", wspeed, sep = ""))
  # dev.off()
  # 
  print(i)
}

outmat <- cbind(tempvec, medvec, maxvec)
colnames(outmat) <- c("mean temp", "median temp", "maximum temp")

write.table(outmat, "outmat-october.txt")


###Parallel approach



getimgs <- function(i)
{
  library(EBImage)
  library(lattice)
  library(viridisLite)
  library(RColorBrewer)
  library(raster)
  
  f <- flist[i]
  
  img <- readImage(f)
  rband <- img[,,1]
  gband <- img[,,2]
  bband <- img[,,3]
  
  rvec <- as.numeric(rband)
  bvec <- as.numeric(bband)
  gvec <- as.numeric(gband)
  
  svec <- rvec + gvec
  ratband <- rband / gband
  
  rdat <- imageData(rband)
  bdat <- imageData(bband)
  gdat <- imageData(gband)
  
  rraster <- raster(rdat)
  braster <- raster(rdat)
  graster <- raster(gdat)
  
  
  ratdat <- imageData(ratband)
  ratraster <- raster(rraster)
  ratvec <- as.numeric(ratdat)
  ratvec[is.infinite(ratvec)] <- .001
  ratvec[is.na(ratvec)] <- .001
  
  #ratvec <- na.exclude(ratvec)
  
  #transvec <- 55.621*ratvec**2 - 590.08*ratvec + 2193.2
  transvec <- 56.552*ratvec**2 -  584.38*ratvec +  2177.7
  transvec[svec < 0.0001] <- 0
  maxval <- 56.552 -  584.38 +  2177.7
  nnval <- quantile(transvec, 0.99)
  
  transvec[transvec > maxval] <- NA
  transsub <- transvec[transvec > 250]
  #transvec[bvec < .01 ] <- 0
  transvec[svec < .075 ] <- 0
  # transvec[transvec <1700] <- 0
  transvec[transvec <250] <- 0
  
  
  transmat <- matrix(transvec, nrow = dim(ratdat)[1], byrow = F)
  #transmat <- t(transmat)
  #transmat <- rotate(transmat)
  #transmat <- rotate(transmat)
  #transmat <- rotate(transmat)
  
  #transmat <- apply(transmat, 1, "rev")
  #transmat <- apply(transmat, 2, "rev")
  
  transimg <- Image(transmat)
  # plot(transimg)
  transraster <- raster(transmat)
  
  transsub <- transvec[transvec > 0]
  meanval <- mean(na.exclude(transsub))
  medval <- median(na.exclude(transsub))
  sdval <- sd(na.exclude(transsub))
  madval <- mad(na.exclude(transsub))
  nval <- length(na.exclude(transsub))
  subraster <- transraster
  subraster[subraster < 500] <- 0
  subraster <- t(subraster)
  transraster <- transraster[transraster]
  fout <- gsub(".tiff", "-contour.png", f)
  png(fout)
  plot(subraster, col = coul, zlim = c(900, maxval))
  dev.off()
  
  outvec <- c(meanval, medval, sdval, madval, nval)
  return(outvec)
}

library(foreach)
#Parallel
require(foreach)
require(parallel)
require(doParallel)

cores <- 6  #32
cl<-makeCluster(cores) #register cores
registerDoParallel(cl, cores = cores)


outlist <- foreach(k=1:nfiles, .combine=rbind, .packages = c("EBImage", "lattice", "viridisLite", "RColorBrewer", "raster")) %dopar% {
  getimgs(k)
}

stopCluster(cl)
closeAllConnections()

plot(outlist[,5])

outsub <- outlist[outlist[,5] > 100,]

plot(outsub[,1])
points(1:nfiles, outlist[,2], col = "red")

nsub <- dim(outsub)[1]

plot(outsub[,1])
points(1:nsub, outsub[,2], col = "red")

plot(outsub[,4])


outmat <- cbind(tempvec, medvec)
colnames(outmat) <- c("mean temp", "median temp")

write.table(outmat, "outmat.txt")
