#########################
# This code sets up and runs the pintail migration network metrics
# Pintail data can be found in metric_inputs_adult_females.xlsx, metric_inputs_adult_males.xlsx, metric_inputs_juvenile_females.xlsx, and metric_inputs_juvenile_males.xlsx
# GR.R does the full G-metric calculation and averaging.


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

seasons <- 3 # Number of seasons or steps in one annual cycle. 
# This must match number of spreadsheets in input files

num_nodes <- 5 # Number of nodes in the network
# This must match the number of initial conditions given in input files

NODENAMES <- c("AK", "PR", "NU", "CA", "GC")
# Used to name row values in CR outputs, should be ordered to match node order in the .xlsx file

NETNAME <- c("adult_female", "adult_male", "juvenile_female", "juvenile_male")  # Give a distinct name for each class as used in input files
# Input files should be in the same directory as RunSpecies.R and named: metric_inputs_NETNAME[[]].xlsx
# Order is important when looking at the code: here we would index [[1]] = class 1 and [[2]] = class 2...

delta <- -1 #Perturbation amount used in G-metric calculation delta=-1 in complete node removal.

PRINT_RESULTS <- TRUE # If true final G-Metric results will print to the screen for the given delta

RUN_PERT <- TRUE # If true the code will rerun for each delta value given in delta_range
# Then the plots G vs delta and GPath vs delta will be created in the directory containing the sourced code -- WARNING this could take some time to run
### RUN_PERT PARAMETERS ###
delta_range <- matrix(c(.5, .4, .3, .2, .1, 0, -.1, -.2, -.3, -.4, -.5, -.6, -.7, -.8, -.9, -1))
node_color <- c("brown", "red", "orange", "blue", "cyan") # AK, PR, NU, CA, GC
NAME <- "Pintails" #Name for plots
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