#' Integrates targets in real data
#'
#' This function integrates the n.ions for each target and calculates Area, Intensity and Extracted Ions Chromatograms for each Target
#' in real data. They are stores in slot R.
#'
#' @export integrateTargets
#' @importFrom tcltk tk_choose.dir
#' @param n.ions number of ions we want to quantify for each targe
#' @examples
#' iTERPSet <- integrateTargets(object=iTERPSet, n.ions = 3)



integrateTargets <- function(object, n.ions=3){
            path <- tk_choose.dir(caption = "Select directory of mzML files to quantify")
            list.of.samples <- list.files(path, pattern=".mzXML", recursive=F)
            object@Sample.name <- list.of.samples
            TIC <- NULL
            EICs <- NULL
            R <- NULL
            for(Smp in 1:length(list.of.samples)){
              cat(Smp)
              sampleRD <- erah:::load.file(paste(path,list.of.samples[Smp], sep="/"))
              mz.val <- c(sampleRD@min.mz:sampleRD@max.mz)
              mz.val.target.spectra <- c(1:600)
              rt <- (1:dim(sampleRD@data)[1])/(sampleRD@scans.per.second*60) + sampleRD@start.time/60
              n.target <- ncol(object@Targets.spectra)
              int <-rowSums(sampleRD@data)
              TIC[[list.of.samples[Smp]]] <- cbind(rt,int)
              for (i in 1:n.target){
                if(is.na(object@Targets.ion.input[[i]])){
                 
                  ###We should ensure that target spectra are within the same range than measured masses
                  ix.target.spectra <- which(as.numeric(rownames(object@Targets.spectra))>=min(mz.val) & as.numeric(rownames(object@Targets.spectra))<=max(mz.val))
                  masses_targets <- mz.val.target.spectra[ix.target.spectra]
                  mass2eics1 <- sort(object@Targets.spectra[ix.target.spectra,i],decreasing = T,index.return=T)$ix[1:n.ions]
                  mass2eics <- masses_targets[mass2eics1]

                  ix.mat <- sapply(mass2eics, function(x) which(mz.val==x))
                  scan.min <- object@Target.comp.scan[i]-object@Targets.scan.window
                  scan.max <- object@Target.comp.scan[i]+object@Targets.scan.window
                  eics <- sampleRD@data[c(scan.min:scan.max),ix.mat]
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




