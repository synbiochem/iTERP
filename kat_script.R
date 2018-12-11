#This is a script for Kat to include new targets in iTERP

# 1-Install development version from SYNBIOCHEM GitHub
if (!("devtools" %in% installed.packages())){
  install.packages("devtools")
}
devtools::install_github("https://github.com/synbiochem/iTERP.git")

# 2-Load iTERP
library(iTERP)

# 3-Load targets
iTERPSet <- NIST2targets(rt.sel = T, ion.sel = T)

# optional step: you might want to plot target spectra for each target
plotMSLibraryTargets(object=iTERPSet,n.target=4, mz.min = 35, mz.max = 200)

# 4- Run all (extract ion chromatograms for selective ions and integrate)
iTERPSet <- runAll(object = iTERPSet, rt.window = 1)

# 5- Retrieve areas
R <- data.frame(t(as.data.frame(lapply(iTERPSet@As, function(x) x))))

