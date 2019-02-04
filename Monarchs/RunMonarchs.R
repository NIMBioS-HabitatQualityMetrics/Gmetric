#########################
# This code sets up and runs the monarch migration network G-metrics
# Monarch data can be found in metric_inputs_monarch.xlsx
# GR.R does the full G calculation and averaging.
# The user must specify which delta they choose, where -1<delta and delta = -1 is complete node removal


## USERS MUST SUPPLY INPUT FILES (metric_inputs_name.xlsx) for each class
## The format of the input files is important
## One tab (sheet) for each season
## Node Attributes at the top to include: Population, Survival Rate, Reproduction or Transition rates, and a list of allowed class transitions
## Then Path Survival rates in matrix form
## Finally Path Transition rates in matrix form

# Clear Workspace
rm(list=ls())

#########################
### USER DEFINED DATA ###
#########################

seasons <- 7 # Number of seasons or steps in one annual cycle. 
# This must match number of spreadsheets in input files

num_nodes <- 4 # Number of nodes in the network
# This must match the number of initial conditions given in input files

NODENAMES <- c("Mexico", "South", "Central", "North")
# Used to name row values in CR outputs, should be ordered to match node order in the .xlsx file

NETNAME <- c("monarch") # Give a distinct name for each class as used in input files
# Input files should be in the same directory as RunSpecies.R and named: metric_inputs_NETNAME[[]].xlsx
# Order is important when looking at the code: here we would index [[1]] = class 1 and [[2]] = class 2

delta <- -1 #Perturbation amount used in G-metric calculations printed to screen delta = -1 means node removal

PRINT_RESULTS <- TRUE # If true final G results will print to the screen

RUN_PERT <- TRUE # If true the code will rerun for each delta value given in delta_range
# Plots G vs delta and GP vs delta will be created in the directory containing RunMonarchs. -- this could take some time to run
### RUN_PERT PARAMETERS ###
delta_range <- matrix(c(.5, .4, .3, .2, .1, 0, -.1, -.2, -.3, -.4, -.5, -.6, -.7, -.8, -.9, -1))
node_color <- c("red", "orange", "blue","green") # MEXICO SOUTH CENTRAL NORTH
NAME <- "Monarchs" # Name for the plots
##########################


#############################################################
### Users should not need to interact with the code below ###
#############################################################


# Set the working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
source("../GR.R")

if(PRINT_RESULTS == TRUE){
  ## Print data to screen
  print("G-metric for all seasons and nodes:")
  print(round(GR,4))
  print("-------------------------------------------------------------------")
  
  print("Season averaged G-metric for all nodes:")
  print(round(GRs,4))
  print("-------------------------------------------------------------------")
  
  print("Network Growth Rate for all seasons:")
  print(round(LAMBDAt,4))
  print("-------------------------------------------------------------------")
  
  print("G-path-metric for all pathways and seasons:")
  print(round(GRP,4))
  print("-------------------------------------------------------------------")
  
  print("Season averaged G-path-metric for all pathways:")
  print(round(GRPs,4))
  print("-------------------------------------------------------------------")
  
}

if(RUN_PERT == TRUE){
  source(paste("../GRplots.R",sep=""))
}