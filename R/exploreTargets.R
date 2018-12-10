
#' Explore targets in real data
#'
#' This function is aimed to establish the rt of target terpenoids in real data
#'
#' @export exploreTargets
#' @importFrom tcltk tk_choose.dir
#' @import erah
#' @param iTERPSet S4 structure class ITERP
#' @param  mz.num number of ions we want to determine for each targe
#' @examples
#' exploreTargets(object=iTERPSet, n.target = 2, n.ions=3, rt.lim=c(2.5,2.8))
#'

exploreTargets <- function(object, n.target = 2, n.ions=3, rt.lim=NULL){
  require(erah)
  old.par <- par(no.readonly = T)
  nom <- tk_choose.files(default = "/Users/mibssmv3/Documents/MARIONA/SYNBIO_PROJECTS/TERPENOIDS_SCREENING",
                         caption = "Select mzXML sample file to explore targets",
                         multi = FALSE);
  titol <- unlist(strsplit(nom, split="/"))[length(unlist(strsplit(nom, split="/")))]
  sampleRD <- erah:::load.file(nom);
  mass2eics <- sort(object@Targets.spectra[,n.target],decreasing = T,index.return=T)$ix[1:n.ions]
  mz.val <- c(sampleRD@min.mz:sampleRD@max.mz)
  ix.mat <- sapply(mass2eics, function(x) which(mz.val==x))
  EICS <- sampleRD@data[,ix.mat]
  rt <- (1:dim(sampleRD@data)[1])/(sampleRD@scans.per.second*60) + sampleRD@start.time/60
  TIC <- rowSums(sampleRD@data)

  if (is.null(rt.lim)==FALSE) {
    ix.rt <- which(rt>rt.lim[1] & rt<rt.lim[2])
    rt <- rt[ix.rt]
    TIC <- TIC[ix.rt]
    ions <- sampleRD@data[,ix.mat]
    EICS <- EICS[ix.rt,]
  }

  par(mfrow=c(1,2))
  plot(rt,TIC, type="l", xlab="rt",ylab="intensity",
  main=paste("TIC-",titol), xlim=rt.lim)
  abline(v =  object@Targets.rt[n.target], col=2,lty=2)
  matplot(rt,EICS, type="l", xlab = "rt", ylab= "intensity",
          main = paste("EICs-",names(object@Targets.spectra)[n.target]), xlim=rt.lim)
  legend("topright", legend = mass2eics, pch=20,
        cex=0.7, col = c(1:length(mass2eics)),bty="n", text.font=10,
        x.intersp=0.1, xjust=1)
  par(old.par)
}
