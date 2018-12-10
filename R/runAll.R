#' Locate Targets in real data
#'
#' This function is aimed to establish the rt of target terpenoids in real data
#' @export runAll
#' @import erah
#' @importFrom tcltk tk_choose.files
#' @param iTERPSet S4 structure class ITERP
#' @param rt.window is the rt window in seconds for maximum width of target peaks
#' @param n.ions are the number of top-intense ions that you want to integrate
#' @examples
#' iTERPSet <- runAll(iTERPSet,rt.window = 3, n.ions=3)
#'

runAll <- function(object, rt.window=3, n.ions=3){
  
  path <- tk_choose.dir(caption = "Select directory of mzML files to quantify")
  list.of.samples <- list.files(path, pattern= ".mzXML", recursive=F)
  list.of.mzXML <- list.files(path, pattern= ".mzXML", recursive=F, full.names = T)
  object@Sample.name <- list.of.samples
  EICs <- NULL
  TICs <- NULL
  As <- NULL
  Is <- NULL
  for (Smp in 1:length(object@Sample.name)){
    cat(Smp)
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
    rt <- hd$retentionTime #retention time in seconds
    tic <- hd$totIonCurrent
    mzR::close(MSfile.open)##close file
    
    eics <- NULL ###EICs list for each sample
    a <- NULL ###Areas list for each sample
    i <- NULL ###Intensities list for each sample
    
    for (target in 1:length(object@Targets.spectra)){
      rt.target <- as.numeric(object@Targets.rt.input[[target]])
      ions.target <- as.numeric(object@Targets.ion.input[[target]])
      if (!is.na(rt.target)){#object contains specified rt
        idx_rt <- which(rt<rt.target*60+rt.window  & rt>rt.target*60-rt.window)
        rt.vector.target <- rt[idx_rt]/60
        if(!is.na(ions.target)){##object contains specified ion
          idx_ion <- which(mz.val==ions.target)
          eics[[target]] <- cbind.data.frame(rt.vector.target, raw.data[idx_rt, idx_ion])
          colnames(eics[[target]]) <- c("rt", "Int")
          a[[target]] <- apply(eics[[target]],2, sum)[2]
          i[[target]] <- apply(eics[[target]],2,max)[2]
          #plot(y=EICS[[3]][[1]], x=rt.vector.target, type="l")
        }else{###object do not contains ions information, we retrieve the n. most intense ions ions EICs in the spectra
          mass2eics1 <- sort(object@Targets.spectra[[target]], decreasing = T, index.return=T)$ix[1:n.ions]
          idx_ion <- sapply(mass2eics1, function(x) which(mz.val==x))
          eics.ion <- lapply(idx_ion, function(ion) {
            eics.ion <- cbind.data.frame(rt.vector.target, raw.data[idx_rt, ion])
            colnames(eics.ion) <- c("rt", "Int")
            return(eics.ion)
          })
          names(eics.ion) <- mass2eics1
          eics[[target]] <- eics.ion
          a.ion <- lapply(eics[[target]], function(x) {apply(x,2, sum)[2]})
          names(a.ion) <- mass2eics1
          i.ion <- lapply(eics[[target]], function(x) {apply(x,2, max)[2]})
          names(i.ion) <- mass2eics1
          a[[target]] <- a.ion
          i[[target]] <- i.ion
        }
      } else {###object do not contains specified rt, we try to locate targets by finding
        #the maximum correlation of target spectra with each scan 
        #Find where in the chromatogram the correlation with spectra is maximized
        spTarg.s <- object@Targets.spectra[[target]][mz.val]
        cor.vectors <- cor(t(raw.data), spTarg.s)
        rt.pred <- rt[which.max(cor.vectors)]##predicted retention time in seconds
        idx_rt <- which(rt<rt.pred+rt.window  & rt>rt.pred-rt.window)
        rt.vector.target <- rt[idx_rt]/60
        object@Targets.rt[target] <- rt.pred/60 ###predicted retention time in minutes
        if(!is.na(ions.target)){##object contains specified ion
          idx_ion <- which(mz.val==ions.target)
          eics[[target]] <- cbind.data.frame(rt.vector.target, raw.data[idx_rt, idx_ion])
          colnames(eics[[target]]) <- c("rt", "Int")
          a[[target]] <- apply(eics[[target]],2, sum)[2]
          i[[target]] <- apply(eics[[target]],2,max)[2]
          #plot(y=EICS[[3]][[1]], x=rt.vector.target, type="l")
        }else{###object do not contains ions information, we retrieve the n. most intense ions ions EICs in the spectra
          mass2eics1 <- sort(object@Targets.spectra[[target]], decreasing = T, index.return=T)$ix[1:n.ions]
          idx_ion <- sapply(mass2eics1, function(x) which(mz.val==x))
          eics.ion <- lapply(idx_ion, function(ion) {
            eics.ion <- cbind.data.frame(rt.vector.target, raw.data[idx_rt, ion])
            colnames(eics.ion) <- c("rt", "Int")
            return(eics.ion)
          })
          names(eics.ion) <- mass2eics1
          eics[[target]] <- eics.ion
          a.ion <- lapply(eics[[target]], function(x) {apply(x,2, sum)[2]})
          names(a.ion) <- mass2eics1
          i.ion <- lapply(eics[[target]], function(x) {apply(x,2, max)[2]})
          names(i.ion) <- mass2eics1
          a[[target]] <- a.ion
          i[[target]] <- i.ion
        }
        
      }
      
    }
    names(eics) <- names(object@Targets.spectra)
    names(a) <- names(object@Targets.spectra)
    names(i) <- names(object@Targets.spectra)
    EICs[[Smp]] <- eics
    As[[Smp]] <- a
    Is[[Smp]] <- i
    TICs[[Smp]] <- cbind.data.frame(rt, tic)
    colnames(TICs[[Smp]]) <- c("rt", "Int")
  }
  names(TICs) <- object@Sample.name
  names(EICs) <- object@Sample.name
  names(As) <- object@Sample.name
  names(Is) <- object@Sample.name
  object@TICs <- TICs
  object@EICs <- EICs
  object@As <- As
  object@Is <- Is
  return(object)
}
            
