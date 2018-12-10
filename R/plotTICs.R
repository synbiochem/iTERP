#' Locate Targets in real data
#'
#' This function is aimed to establish the rt of target terpenoids in real data
#' @export plotTICs
#' @param iTERPSet S4 structure class ITERP
#' @param sample2plot is  the number of dample to plot. If it is null an overlay of all TICS is iTERPset
#' is going to be displayed
#' @param rt.lim is the rt.limit you want to plot
#' @examples
#' plotTICs(object=iTERPSet, sample2plot=6,rt.lim=c(2.5, 2.8))
#' plotTICs(object=iTERPSet,rt.lim=c(2.5, 2.7), sample.class=NULL)


plotTICs <- function(object, sample2plot=NULL,  rt.lim=NULL, sample.class=NULL){
  old.par <- par(no.readonly = T)
  minTICs <- min(sapply(object@TIC, function(x) dim(x))[1,])
  t <- do.call(cbind, lapply(object@TIC,function(x) x[c(1:minTICs),"rt"]))
  int <- do.call(cbind, lapply(object@TIC,function(x) x[c(1:minTICs),"int"]))

  if(!is.null(sample2plot)==TRUE){
    if(!is.null(rt.lim)==TRUE){
      int.s <- int[,sample2plot]
      t.s <- t[,sample2plot]
      max.int <- max(int.s[t.s>rt.lim[1] & t.s<rt.lim[2]])
      int.lim <- c(0,max.int)
    }else{
      int.lim <- NULL
    }
    plot(t[,sample2plot],int[,sample2plot], type="l",
         xlab="rt", ylab="counts", xlim=rt.lim, ylim=int.lim,
         main=paste("TIC-",object@Sample.name[sample2plot],sep="")
    )
  }else{

    if(!is.null(rt.lim)==TRUE){
      max.int <- max(int[t>rt.lim[1] & t<rt.lim[2]])
      int.lim <- c(0,max.int)
    }else{
      int.lim <- NULL
    }
    matplot(t,int, type="l", xlab = "rt",
            ylab= "intensity",
            main = "Overlayed TIC", xlim=rt.lim, ylim=int.lim,
            col=sample.class)
  }
  par(old.par)###preserves parameters
}





