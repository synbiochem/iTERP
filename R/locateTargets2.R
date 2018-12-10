
#' Locate Targets in real data
#'
#' This function is aimed to establish the rt of target terpenoids in real data
#' @export locateTargets
#' @import erah
#' @importFrom tcltk tk_choose.files
#' @param iTERPSet S4 structure class ITERP
#' @param scan.window is the rt window in scans for maximum width of target peaks
#' @examples
#' iTERPSet <- locateTargets(iTERPSet,scan.window = 10)
#'

locateTargets <- function(object, scan.window=10){
            #rt.window is the scan window + - where peaks are located

            nom <- tk_choose.files(default = "/Users/mibssmv3/Documents/MARIONA/SYNBIO_PROJECTS/TERPENOIDS_SCREENING",
                                   caption = "Select mzXML sample file to locate targets",
                                   multi = FALSE);
            titol <- unlist(strsplit(nom, split="/"))[length(unlist(strsplit(nom, split="/")))]
            require(erah)
            sampleRD <- erah:::load.file(nom);
            object@Sample.to.locate.targets <- titol;
            spTarg.s <- object@Targets.spectra[sampleRD@min.mz:sampleRD@max.mz,]
            ###Find rt
            no.scans <- dim(sampleRD@data)[1]
            rt <- (1:no.scans)/(sampleRD@scans.per.second*60) + sampleRD@start.time/60
            
            rt.mean <- NULL
            rt.vector1 <-NULL
            rt.vector2 <-NULL
            CmpScan2 <- NULL
            rt.met2 <- NULL
            #Find where in the chromatogram the correlation with spectra is maximized
            cor.vectors <- cor(t(sampleRD@data), spTarg.s)
            colnames(cor.vectors) <- colnames(spTarg.s)
            
            ####Locate rt of input targets in the sample
            for (target in 1:length(object@Targets.rt.input)){
              rt.target <- as.numeric(object@Targets.rt.input[[target]])
              if (!is.na(rt.target)){#Tenemos rt como input
                 CmpScan <- which(abs(rt-rt.target)==min(abs(rt-rt.target)))
                 mat <- sampleRD@data[(CmpScan-scan.window):(CmpScan+scan.window),]
                 max_tic <- which.max(rowSums(mat));#Put the peak centre where TIC is maximum
                 rt.vector1[[target]] <- rt[(CmpScan-scan.window):(CmpScan+scan.window)]
                 rt.vector2[target] <- rt.vector1[[target]][max_tic]
                 rt.met2[[target]]<- rt.vector1
                 CmpScan2[target] <- CmpScan
              }else{
                Kstart <- which.max(cor.vectors[,target]) - scan.window
                Kend <- which.max(cor.vectors[,target]) + scan.window
                rt.vector1[[target]] <- rt[Kstart:Kend]
                rt.mean[target] <- which.max(cor.vectors[,target])/(sampleRD@scans.per.second*60) + sampleRD@start.time/60
                CmpScan <- (rt.mean[target] - sampleRD@start.time/60)*(sampleRD@scans.per.second*60)
                mat <- sampleRD@data[(CmpScan-scan.window):(CmpScan+scan.window),]
                max_tic <- which.max(rowSums(mat));#Put the peak centre where TIC is maximum
                rt.vector2[target] <- rt.vector1[[target]][max_tic]
                CmpScan2[target] <- (rt.vector2[target] - sampleRD@start.time/60)*(sampleRD@scans.per.second*60)
                Kstart2 <- CmpScan2[target] - scan.window
                Kend2 <- CmpScan2[target] + scan.window
                rt.met2[[target]] <- rt[Kstart2:Kend2]
                }
            }

            Targets.rt <- rt.vector2; names(Targets.rt) <- colnames(object@Targets.spectra)
            Targets.rt.vectors <- rt.met2; names(Targets.rt.vectors) <- colnames(object@Targets.spectra)
            object@Targets.rt <- Targets.rt
            object@Targets.rt.vectors <- Targets.rt.vectors
            object@Target.comp.scan <- CmpScan2
            object@Targets.scan.window <- scan.window
            #Plot TIC together with central peaks found
            TIC_i <- cbind(((1:no.scans)/(sampleRD@scans.per.second*60) + sampleRD@start.time/60),rowSums(sampleRD@data))
            plot(x=TIC_i[,1],y=TIC_i[,2], type="l", xlab="rt", ylab="counts",
                 main=paste("TIC-", titol, sep=""))
            abline(v =  object@Targets.rt, col=2,lty=2)
            return(object)
          }
