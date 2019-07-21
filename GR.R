########## G-METRIC CALCULATION CODE ############
# This code calculates:
# GR - the G-metric for each class and each node at each season step
# GRt - the class population weighted G-metric for each node at each season step
# GRs - the season population weighted G-metric for each node across the annual cycle
# LAMBDAt - the network growth rate for each season
# GAMMAt - the network growth rate for each season after applying the delta impact (removal, degradation, or enhancement)

library(XLConnect)

### Users should not need to interact with the code below ###

AMATRIX <- list()
AMATRIX_CR <- list()
f_update <- list()
p_update <- list()
s_update <- list()
rates <- c("population","survival","reproduction","transition")

### Find the number of classes
NUMNET <- length(NETNAME)

# Set the input file names
input_file_names <- matrix(0,1,NUMNET)
for (i in 1:NUMNET){
  input_file_names[i]<-c(paste("metric_inputs_",NETNAME[i],".xlsx", sep="")) 
}

#####################################
### CHECK VALIDITY OF INPUT FILES ###
#####################################

for (i in 1:NUMNET){
  CHECK_WORKBOOK <- file.exists(input_file_names[i])
  if (CHECK_WORKBOOK==F){stop(paste("The Workbook,", input_file_names[i],"could not be found. \n *** Please check that the end of the input file names match NETNAME"))}
}
rm(i)
for (i in 1:NUMNET){
  SPREADSHEETS <- getSheets(loadWorkbook(input_file_names[i]))
  NUMSHEETS <- length(SPREADSHEETS)
  if (NUMSHEETS != seasons){stop(paste("\n The Workbook,", input_file_names[i],"has an incorrect number of sheets. \n *** It should contain one sheet for initial conditions and one sheet for each season."))}
}

###################################
### IMPORT NODE CHARACTERISTICS ###
###################################

NODE_ATTRIBUTES <- list()
PATH_ATTRIBUTES <- list()
TRANS_ATTRIBUTES <- list()

## Node data belongs on the top starting in cell C3
## Path data should be contained in nXn ranges below the node data with 3 rows separating
## Survival probabilities should start in cell C(n+6)
## Transition rates should start in cell C(2n+9) where n=number of nodes


for (i in 1:NUMNET) {
  NODE_ATTRIBUTES[[NETNAME[i]]]<- list()
  PATH_ATTRIBUTES[[NETNAME[i]]]<- list()
  TRANS_ATTRIBUTES[[NETNAME[i]]] <- list()
  for (k in 1:seasons){
    # READ IN THE DATA
    SR <- 2 # row in the spreadsheet where we start reading node data
    SC <- 3 # column in the spreadsheet where we start reading node data
  
    ###____seasonal node characteristics ____###
    NODE_ATTRIBUTES[[NETNAME[i]]][[k]] <- list()
    for (m in 1:3) {
      NODE_ATTRIBUTES[[NETNAME[i]]][[k]][[rates[[m]]]]<- matrix(0,nrow = num_nodes,ncol = 1)
      NODE_ATTRIBUTES[[NETNAME[i]]][[k]][[rates[[m]]]]<- readWorksheetFromFile(input_file_names[i], sheet=k, startRow=SR, endRow=SR+num_nodes, startCol=SC, endCol=SC)
      SC <- SC + 1
      # Check that there are no empty cells in node characteristics
      if(any(is.na(NODE_ATTRIBUTES[[NETNAME[i]]][[k]][[m]]))){
        stop(paste('\n Check the input file! \n *** There is at least one empty node characteristic in season', k))}
    }
    
    TRANS_ATTRIBUTES[[NETNAME[i]]][[k]]<- matrix(0,nrow=NUMNET,ncol=1)
    TRANS_ATTRIBUTES[[NETNAME[i]]][[k]]<- readWorksheetFromFile(input_file_names[i], sheet=k, startRow=SR, endRow=SR+NUMNET, startCol=7, endCol=7)
    # Check that the transition atributes are in the correct order
    temp <- matrix(0,nrow=NUMNET,ncol=1)
    temp <- readWorksheetFromFile(input_file_names[i], sheet=k, startRow=SR, endRow=SR+NUMNET, startCol=6, endCol=6)
    if(any(temp != NETNAME)){
      stop(paste('\n Check the input file! \n *** The order of classes for allowed transitions must match the network order in season', k, 'for class', NETNAME[i],'\n *** check the variables temp and NETNAME.'))
    }
    
    SC <- 3
    ###____seasonal path characteristics ____###
    PATH_ATTRIBUTES[[NETNAME[i]]][[k]]<- list()
    SR <- num_nodes + 5
    for (n in 1:2){
      PATH_ATTRIBUTES[[NETNAME[i]]][[k]][[n]] <- readWorksheetFromFile(input_file_names[i],sheet=k,startRow=SR,endRow=SR+num_nodes,startCol=SC,endCol= SC + num_nodes)
      SR <- (2 * num_nodes) + 8
    }
  } 
}

#####################
### CREATE N LIST ###
#####################
N <- matrix(0, nrow = NUMNET*num_nodes, ncol = seasons)
for (SN in 1:seasons) {
  for (ND in 1:num_nodes) {
    for (CL in 1:NUMNET) {
      N[CL+NUMNET*(ND-1),SN]<- NODE_ATTRIBUTES[[NETNAME[[CL]]]][[SN]]$population[[1]][ND]
    }
  }
}

#########################
### CREATE A MATRICES ###
#########################

for (SN in 1:seasons) {
## MODEL FUNCTIONS ##
  f_update[[SN]]<- list()
  for (i in 1:num_nodes){
    node <- i
    f_temp <- matrix(0,nrow = NUMNET,ncol = NUMNET)
    for (j in 1:NUMNET) {
      for (k in 1:NUMNET) {
        if (j==k){
          # Choose from survival rates: NODE_ATRIBUTES[[NETNAME]][[season]][[2=survival]][[1]][[node=i]] 
          f_temp[j,k]<- NODE_ATTRIBUTES[[NETNAME[[j]]]][[SN]][[2]][[1]][[i]]
        }
        else{
          # Choose from transition/reproduction rates and then test if transitions are allowed (0=no, 1=yes)
          f_temp[j,k]<- NODE_ATTRIBUTES[[NETNAME[[j]]]][[SN]][[3]][[1]][[i]]*TRANS_ATTRIBUTES[[NETNAME[[j]]]][[SN]][[1]][[k]]
        }
     }
    }
    f_update[[SN]][[i]] <- f_temp
  }
  rm(i)
  p_update[[SN]] <- list()
  s_update[[SN]] <- list()
  for (i in 1:NUMNET){
    type <- i
    p_update[[SN]][[i]] <- PATH_ATTRIBUTES[[i]][[SN]][[2]]   
    s_update[[SN]][[i]] <- PATH_ATTRIBUTES[[i]][[SN]][[1]]
  }
rm(i)


## MATRICEES NEEDED FOR UPDATE ##
## F Block Diagonal ##
FBLOCK <- matrix(0,nrow = NUMNET*num_nodes, ncol = NUMNET*num_nodes)
for (i in 1:num_nodes){
  E <- matrix(0,nrow = num_nodes, ncol = num_nodes)
  E[i,i] <- 1
  FBLOCK <- FBLOCK + kronecker(E, t(matrix(f_update[[SN]][[i]],nrow=NUMNET,ncol=NUMNET)))
}
rm(i,E)
## Q Block Diagonal ##
QBLOCK <- matrix(0,nrow = NUMNET*num_nodes, ncol = NUMNET*num_nodes)
for (i in 1:NUMNET){
  E <- matrix(0,nrow = NUMNET,ncol = NUMNET)
  E[i,i] <- 1
  QBLOCK <- QBLOCK + kronecker(t(matrix(unlist(s_update[[SN]][[i]]),nrow = num_nodes,ncol = num_nodes) * matrix(unlist(p_update[[SN]][[i]]),nrow = num_nodes,ncol = num_nodes)),E)
}
rm(i,E)

AMATRIX <- QBLOCK %*% FBLOCK
AMATRIX_CR[[SN]] <- AMATRIX
AMATRIX_CR[[SN+seasons]] <- AMATRIX
}

################################
### CALCULATE GR NODE VALUES ###
################################

# First calculate per-capita contribution matrix CR - used to find network growth rate
POPnode <- matrix(0,nrow = num_nodes, ncol = seasons)
CR <- matrix(0,nrow = nrow(AMATRIX), ncol = seasons) #Cr values for each class/node/season
ONES <- matrix(1,nrow = nrow(AMATRIX), ncol = 1)
#Choose each season as a focal season
for (SN in 1:seasons){
  CRtemp <- diag(nrow(AMATRIX))
  taustart <- SN  
  tauend <- SN+seasons-1
  for (i in taustart:tauend){
    CRtemp <- CRtemp %*% t(matrix(AMATRIX_CR[[i]], nrow = nrow(AMATRIX)))
  }
  CR[,SN] <- CRtemp %*% ONES
  for (j in 1:num_nodes){
    POPnode[j,SN] <- sum(N[((j-1)*NUMNET+1):(j*NUMNET),SN])
  }
}



## Name columns and rows
SEASONNAMES <- c()
for (i in 1:seasons){
  SEASONNAMES <- c(SEASONNAMES, paste("season", i))
}
NAMESofROW <- c()
num_class <- length(NETNAME)
for (i in 1:num_nodes){
  for (j in 1:num_class){
    NAMESofROW <- c(NAMESofROW, paste(NODENAMES[i], NETNAME[j]))
  }
}
NAMESofROW <- gsub("_","",NAMESofROW)



## Population Proportion and Network Growth Rate
LAMBDAt <- matrix(0, 1, seasons)
WR <- matrix(0, NUMNET*num_nodes, seasons) #Contains population porportion for each class at each node.
WRt <- matrix(0,num_nodes,seasons) #Contains the population proportiona for each node, summing the classes
WRs <- matrix(0,num_nodes,1) #Contains the network wide average WR values - divide by number of seasons
for (SN in 1:seasons){
  TIME <- SN
  WR[,SN] <- N[,TIME]/sum(N[,TIME])
  TIME <- TIME + 1
}
colnames(WR) <- SEASONNAMES
rownames(WR) <- NAMESofROW

for (i in 1:num_nodes){
  if(NUMNET==1){
    WRt[i,] <- WR[i,]
  }
  else{
    WRt[i,] <- colSums(WR[((i-1)*NUMNET+1):(i*NUMNET),])
  }
}

WRs <- rowSums(WRt)/3

for (SN in 1:seasons){
  LAMBDAt[SN] <- t(WR[,SN])%*%CR[,SN]
}
rownames(LAMBDAt)<-c("network growth rate")
colnames(LAMBDAt)<-SEASONNAMES


## G-metric  and network growthrate in absense of node r

GR <- matrix(0,num_nodes,seasons)
GRs <- matrix(0,num_nodes,1) # Season averaged G-metric
GAMMAt <- matrix(0,num_nodes,seasons)

Inc <- diag(1, nrow(AMATRIX), ncol(AMATRIX))
for (j in 1:num_nodes){
  for (SN in 1:seasons){
    Drr <- matrix(0, nrow(AMATRIX), ncol(AMATRIX))
    Enc <- matrix(0, nrow(AMATRIX), ncol(AMATRIX))
    xstart <- j*NUMNET-NUMNET+1
    xend <- j*NUMNET
    for (x in xstart:xend){
      TEMP <- matrix(0,nrow(AMATRIX), ncol(AMATRIX))
      TEMP[x,x] <- delta
      Enc <- Enc + TEMP  
    }
    Drr <- Inc+Enc
    CRnewtemp <- diag(nrow(AMATRIX))
    taustart <- SN  
    tauend <- SN+seasons-1
    for (i in taustart:tauend){
      CRnewtemp <- CRnewtemp %*% t(matrix(AMATRIX_CR[[i]], nrow = nrow(AMATRIX)))%*%Drr
    }
    GAMMAt[j,SN] <- WR[,SN]%*%Drr%*%CRnewtemp %*% ONES
  }
}

for (j in 1:num_nodes){
  GR[j,] <- abs(LAMBDAt - GAMMAt[j,])  
}
rownames(GR) <- NODENAMES
colnames(GR) <- SEASONNAMES

POPsums <-colSums(POPnode)
TEMP <- matrix(0,num_nodes,seasons)
for (i in 1:seasons){
  TEMP[,i] <- matrix(POPsums[i],1,num_nodes)
}
GRs <- matrix(rowSums(GR*TEMP)/sum(POPsums),nrow = num_nodes,ncol = 1)

rownames(GRs) <- NODENAMES
colnames(GRs) <- c("G-metric")


########### Pathway G-metric ###############

GRP <- matrix(0,num_nodes*num_nodes,seasons)
GRPs <- matrix(0,num_nodes*num_nodes,1)
GAMMAPt <- matrix(0,num_nodes*num_nodes,seasons)

Ic <- matrix(1, nrow=NUMNET, ncol=NUMNET)
In <- matrix(1, nrow=num_nodes, ncol= num_nodes)
count <- 1

for (i in 1:num_nodes){
  for (j in 1:num_nodes){
    if(i+j == 2){
      PATHNAMES <- c("$G_{11,t}$")
    }else{
      PATHNAMES <- c(PATHNAMES,paste("$G_{",i,j,",t}$",sep=""))
    }
    
    # Consider each possible connection ij between each of the nodes
    for (SN in 1:seasons){
      # Consider each season for the focal node
      Brd <- matrix(0,  nrow(In), ncol(In))
      Enrd <- matrix(0,  nrow(In), ncol(In))
      
      Enrd[i,j] <- delta
      Brd <- In+Enrd
      
      TEMP <- diag(nrow(AMATRIX)) # Start with an identity matrix ncXnc
      taustart <- SN
      tauend <- SN+seasons-1
      for (time in taustart:tauend){
        TEMP <- TEMP %*% (t(matrix(AMATRIX_CR[[time]], nrow = nrow(AMATRIX))) * kronecker(Brd,Ic))
      }
      GAMMAPt[count,SN] <- WR[,SN]%*%TEMP%*%ONES
    }
    count <- count + 1
  }
}

rownames(GAMMAPt) <- PATHNAMES
colnames(GAMMAPt) <- SEASONNAMES

nc <- num_nodes*num_nodes
for (j in 1:nc){
  for (SN in 1:seasons){
    GRP[j,SN] <- abs(LAMBDAt[SN] - GAMMAPt[j,SN]) # CHANGED THIS to absolute vales to correct for negative criticality
  }
}
rownames(GRP) <- PATHNAMES
colnames(GRP) <- SEASONNAMES

POPsums <-colSums(POPnode)
TEMP <- matrix(0,num_nodes*num_nodes,seasons)
for (i in 1:seasons){
  cols <- num_nodes*num_nodes
  TEMP[,i] <- matrix(POPsums[i],1,cols)
}
GRPs <- matrix(rowSums(GRP*TEMP)/sum(POPsums),num_nodes*num_nodes,1)
rownames(GRPs) <- PATHNAMES
colnames(GRPs) <- c("G-path Season Average")

