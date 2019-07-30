This code calculates the G-metric as defined in the paper:

Quantifying the contribution of habitats and pathways to a spatially structured population facing environmental change. (2019)
by Christine Sample, Joanna A. Bieri, Benjamin Allen, Yulia Dementieva, Alyssa Carson, Connor Higgins, Sadie Piatt, Shirley Qiu, Summer Stafford, Brady J. Mattsson, Darius Semmens, James E. Diffendorfer, and Wayne E. Thogmartin

Code writen, developed, and tested by: Joanna Bieri, Christine Sample, and Summer Stafford.

NEED TO INSTALL R LIBRARIES: XLconnect

To run the code for any of the case studies:

1. Monarchs
2. Pintails
3. Pulliam Metapopulation
4. The hypothetical files are included to give another example of how to set up the network within the G-metric code, but this network and these results are not included in the paper

Navigate to the example folder.

Here you will find a spreadsheet (.xlsx) that contains data for each class in the example. The formating and location of these files is important since the G code reads in the data based on location in the spreadsheet. There should be one tab for each season in the annual cycle. For each season, Node Attributes are at the top of the spreadsheet, followed by Path Transitions and then Path Survival rates.

Also within each folder is the Run<speciesname>.R code. This is the code that should be sourced to get the G results. The User Defined Data in this code should match what is provided in the spreadsheets.

The main mathematical calculation is found in GR.R and code for perturbing the system and plotting results is found in GRplots.R. Users should not need to interact with this code.
