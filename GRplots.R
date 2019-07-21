############## PERTURBATIONS ON DELTA - USES GR.R G-metric CODE ####
# This code reruns the G-metric code for all values stored in delta_range, defined in runspecies.R
# After running and storing the G-metric node and pathway values plots are created of the data.
# G-metric values are stored in GR_VALUES
# Gpathway values are stored in GRP_VALUES
# The plots are created and stored in the working directory as .pdf files.

GR_VALUES <- list()
GRs_VALUES <- list()
GRP_VALUES <- list()
GRPs_VALUES <- list()
num <- 1

PLOTGR <- matrix(0,nrow=num_nodes,ncol=nrow(delta_range))
PLOTGAMMA <- matrix(0,nrow=num_nodes,ncol=nrow(delta_range))
PLOTGRP <- matrix(0,nrow=num_nodes*num_nodes,ncol=nrow(delta_range))
PLOTGAMMAP <- matrix(0,nrow=num_nodes*num_nodes,ncol=nrow(delta_range))

# For each delta values in the range find the G values and save them for plotting
print("Calculating G for a range of delta perturbations:")
for (DEL in delta_range){
  delta <- DEL
  # CALCULATE G for each given delta
  print(paste( "delta = ", DEL))
  source(paste("../GR.R",sep=""))
  
  # Save G values for each delta
  GR_VALUES[[num]] <- GR
  GRs_VALUES[[num]] <- GRs
  GRP_VALUES[[num]] <- GRP
  GRPs_VALUES[[num]] <- GRPs
  
  PLOTGR[,num] <- GRs 
  PLOTGRP[,num] <- GRPs
  num <- num + 1
  
}

# Find gamma where sign(delta_range) corrects for the absolute value
pathnum <- 1
for(i in 1:num_nodes){
  PLOTGAMMA[i,] <- 1+sign(delta_range)*PLOTGR[i,]
  for(j in 1:num_nodes){
    PLOTGAMMAP[pathnum,] <- 1+sign(delta_range)*PLOTGRP[pathnum,]
    pathnum <- pathnum + 1
  }
}

print("Preparing fountain plots.")
### PLOTS ###

# NODES #
#
#


pdf(paste(NAME,"_GvsDelta_plots.pdf", sep="")) ## START PLOT - save as .pdf file
symb <- 0
plot(c(0),c(0),
     type="l", 
     main=(paste("Contribution as a function of perturbation ", delta)),
     ylim = c(min(PLOTGR),max(PLOTGR)), 
     xlim = c(min(delta_range),max(delta_range)), 
     xlab=expression(delta),
     ylab="",
     cex.lab=1.75, 
     cex.axis=1.75, 
     cex.main=1.75
     )
axis(side = 1, at = seq(-1,0.5, by=.1), tck = 1, lty = 2, col = "grey", labels = NA)
axis(side = 2,  tck=1, lty=2, col = "grey", labels = NA)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
par(xpd = TRUE) #Draw outside plot area
text(x = corners[1]-.19, y = mean(corners[3:4]), expression( italic(G[r]) ), cex=1.75)
for(i in 1:num_nodes){
  par(new=TRUE)
  yplot <- PLOTGR[i,]
  plot(delta_range,yplot,
       type="o", 
       lwd=1, 
       pch=symb, 
       col=node_color[symb+1],
       ylim = c(min(PLOTGR),max(PLOTGR)), 
       xlim = c(min(delta_range),max(delta_range)), 
       axes = FALSE,  
       xlab="", 
       ylab="")
  symb <- symb +1
  
}
dev.off() ## END PLOT

# PATHWAYS #
#
#

pdf(paste(NAME,"_GRPathvsDelta_plots.pdf", sep="")) ## START PLOT - save as .pdf file
symb <- 0
plot(c(0),c(0),
     type="l", 
     main=expression(paste("Pathway Contribution, as a function of ", delta)),
     ylim = c(min(PLOTGRP),max(PLOTGRP)), 
     xlim = c(min(delta_range),max(delta_range)), 
     xlab=expression(delta),
     ylab="",
     cex.lab=1.75, 
     cex.axis=1.75, 
     cex.main=1.75)
axis(side = 1, at = seq(-1,0.5, by=.1), tck = 1, lty = 2, col = "grey", labels = NA)
axis(side = 2,  tck=1, lty=2, col = "grey", labels = NA)
corners = par("usr") #Gets the four corners of plot area (x1, x2, y1, y2)
par(xpd = TRUE) #Draw outside plot area
text(x = corners[1]-.19, y = mean(corners[3:4]), expression( italic(G[r]) ), cex=1.75)
end <- num_nodes*num_nodes
for(i in 1:end){
  if(PLOTGRP[i,1]>10^(-5)){
    par(new=TRUE)
    yplot <- PLOTGRP[i,]
    plot(delta_range,yplot,
         type="o", 
         lwd=1, 
         pch=symb, 
         col=node_color[ceiling(i/num_nodes)],
         ylim = c(min(PLOTGRP),max(PLOTGRP)), 
         xlim = c(min(delta_range),max(delta_range)), 
         axes = FALSE,  
         xlab="", 
         ylab="")
    symb <- symb +1
  }
}
dev.off() ## END PLOT
