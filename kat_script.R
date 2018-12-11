#This is a script for Kat to include new targets in iTERP

# 1-Install development version from SYNBIOCHEM GitHub
if (!("devtools" %in% installed.packages())){
  install.packages("devtools")
}
devtools::install_github("https://github.com/synbiochem/iTERP.git")


setwd("/Volumes/shared/Maria/181210_ReferenceFiles_iTERP")

# 2-Load iTERP
library(iTERP)

# 3- I have done the first step (load targets) for you. You might want to load the following 
# workspace containing and iTERPSet object with targets spectral information
load("/Volumes/shared/Maria/181210_ReferenceFiles_iTERP/iTERPSet_new_targets.RData")

#Alternatively you might want to repeat it:  
#iTERPSet <- NIST2targets(rt.sel = T, ion.sel = T) #You will be asked to redirect to the files where .MSP files for pure standards are allocated

# optional step: you might want to plot target spectra for each target
plotMSLibraryTargets(object=iTERPSet,n.target=4, mz.min = 35, mz.max = 200)

# 4- Run all. You will be asked for the directory where mzXML files are. Extract ion chromatograms for selective ions and integration are calaculted for all mzXML files in the directory
iTERPSet <- runAll(object = iTERPSet, rt.window = 1)

# 5- Retrieve areas
R <- data.frame(t(as.data.frame(lapply(iTERPSet@As, function(x) x))))
