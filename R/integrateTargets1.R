
#' Integrates targets in real data
#'
#' This function integrates the n.ions for each target and calculates Area, Intensity and Extracted Ions Chromatograms for each Target
#' in real data. They are stores in slot R.
#'
#' @export integrateTargets1
#' @importFrom tcltk tk_choose.dir
#' @param n.ions number of ions we want to quantify for each targe
#' @examples
#' iTERPSet <- integrateTargets(object=iTERPSet, n.ions = 3)



integrateTargets1 <- function(object, n.ions=3){
            path <- tk_choose.dir(caption = "Select directory of mzML files to quantify")
            list.of.samples <- list.files(path, pattern=".mzXML", recursive=F)
            list.of.mzXML <- list.files(path, pattern=".mzXML", recursive=F, full.names = T)
            object@Sample.name <- list.of.samples
            TIC <- NULL
            EICs <- NULL
            R <- NULL
            for(Smp in 1:length(list.of.samples)){
              cat(Smp, list.of.samples[Smp], sep="")
              require(mzR)
              MSfile.open <- mzR::openMSfile(filename = list.of.mzXML[Smp])
              hd <- mzR::header(MSfile.open)
              scans.profile <- hd$seqNum
              lowMZ <- round(min(hd$lowMZ))
              highMZ <- round(max(hd$highMZ))
              mz.targets <- unique(sort(as.numeric(unlist(lapply(object@Targets.spectra, function(x) {which(x!=0)})))))
              mz.targets.scan <- mz.targets[which(mz.targets>=lowMZ & mz.targets<=highMZ)]
              raw.data <- mzR::get3Dmap(object=MSfile.open, scans=scans.profile,
                                        lowMz=min(mz.targets.scan), highMz=max(mz.targets.scan), resMz = 1)
              mz.val <- c(min(mz.targets.scan):max(mz.targets.scan))
              
              n.target <- ncol(object@Targets.spectra)
              TIC[[list.of.samples[Smp]]] <- hd$totIonCurrent
              mass.targets <- lapply(object@Targets.spectra, function(x) {which(x!=0)})
              for (i in 1:n.target){
                if(is.na(object@Targets.ion.input[[i]])){###ions are not specified therefore untargeted analytics
                  ####we will find automatically the spectra of targets by finding the maximum correlation
                  #with raw data matrix and set up automatically rt for those targets and we will integrate the 
                
                  mass2eics1 <- sort(object@Targets.spectra[[i]], decreasing = T, index.return=T)$ix[1:n.ions]
                  ix.mat <- sapply(mass2eics1, function(x) which(mz.val==x))
                  scan.min <- object@Target.comp.scan[i]-object@Targets.scan.window
                  scan.max <- object@Target.comp.scan[i]+object@Targets.scan.window
                  eics <- raw.data[c(scan.min:scan.max),ix.mat]
                  colnames(eics) <- mass2eics
                  rownames(eics) <- rt[c(scan.min:scan.max)]
                  A <- apply(eics,2,sum)
                  I <- apply(eics,2,max)
                  r <- list(eics,A,I); names(r) <- c("EICS","A","I")
                } else {
                  mass2eics <- object@Targets.ion.input[[i]]
                  ix.mat <- which(mz.val==mass2eics)
                  scan.min <- object@Target.comp.scan[i]-object@Targets.scan.window
                  scan.max <- object@Target.comp.scan[i]+object@Targets.scan.window
                  eics <- as.data.frame(sampleRD@data[c(scan.min:scan.max),ix.mat])
                  colnames(eics) <- mass2eics
                  rownames(eics) <- rt[c(scan.min:scan.max)]
                  A <- apply(eics,2,sum)
                  I <- apply(eics,2,max)
                  r <- list(eics,A,I); names(r) <- c("EICS","A","I")
                } 
                R[[i]] <- r
              }
            names(R) <- colnames(object@Targets.spectra)
            object@R[[Smp]] <- R
            }
            names(object@R) <- list.of.samples
            object@TIC <- TIC
            return(object)
}




