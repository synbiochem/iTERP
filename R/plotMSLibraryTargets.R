#' Locate Targets in real data
#'
#' This function is aimed to establish the rt of target terpenoids in real data
#' @export
#' @param iTERPSet S4 structure class ITERP
#' @param scan.window is the rt window in scans for maximum width of target peaks
#' @examples
#' iTERPSet <- locateTargets(iTERPSet,scan.window = 10)
#'




plotMSLibraryTargets <- function(object,n.target=2, mz.min=35, mz.max=300){
            #n.target=target number to be plot
            #spTarg=targets data with target spectra as colums and rows 1:600 representing mass measured
            #mz.min, mzmax= range to plot
            spTarg <- object@Targets.spectra
            rel.int <- function(x=spec){(x/(max(x)))*100}
            spec <- rel.int(x=spTarg[c(mz.min:mz.max),n.target])
            mz <- c(mz.min:mz.max);
            labels <- rep("",times=length(mz))
            labels[which(spec>10)] <- mz[which(spec>10)]

            plot(mz, spec, type="h",
                 xlab = "m/Z", ylab = "Rel. Intensity",
                 main=colnames(spTarg[n.target]), ylim=c(0, 110))
            text(x = mz, y = spec,
                 label = labels, pos = 3, cex = 0.8, col = "red")
          }
