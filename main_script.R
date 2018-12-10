#######This is to install package#######
devtools::install()

#This is to document functions
devtools::document()

######This is to built binary############
devtools::build(binary = FALSE)

#Built package site
pkgdown::build_site()

####This is to intall the package from gitlab#############
devtools::install_github("https://github.com/synbiochem/iTERP.git")

#######This is to execute the first version of iTERP########
library(iTERP)

iTERPSet <- NIST2targets(rt.sel = T, ion.sel = T)
iTERPSet <- NIST2targets(rt.sel = T, ion.sel = F)
##Plot mass spectra from target compounds in library
#plotMSLibraryTargets(object=iTERPSet,n.target=1, mz.min = 35, mz.max = 160)

#Locate targets in test sample
iTERPSet <- locateTargets(object = iTERPSet,scan.window = 5)

####Explore wheter these targets appear or not
exploreTargets(object=iTERPSet, n.target = 2, n.ions=3, rt.lim=c(2,3.2))

#Integrate targets
iTERPSet <- integrateTargets(object=iTERPSet, n.ions = 3)
plotEICs(object = iTERPSet, n.target = 2)

#Retrieve Data
resultsA <- getResults(object=iTERPSet, type=c("A"))


###Calculating Kovats Index

calKovatsR <- calKovats()
KI <- KovatsIndex(calibration=calKovatsR, rt_unk)


####Runing runAll runAll executes locate targets and integrate targets all in one. No need for selecting a sample locating targets
library(iTERP)
iTERPSet <- NIST2targets(rt.sel = T, ion.sel = T)
iTERPSet <- runAll(object = iTERPSet, rt.window = 1)

target <- 7
rt <- lapply(iTERPSet@EICs, function(x) x[[target]]$rt)
int <- lapply(iTERPSet@EICs, function(x) x[[target]]$Int)
plot(unlist(rt), unlist(int), type="l")


R <- data.frame(t(as.data.frame(lapply(iTERPSet@As, function(x) x))))
