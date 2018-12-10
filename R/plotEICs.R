#' plotEICs
#'
#' This function is aimed to establish the rt of target terpenoids in real data
#' @export plotEICs
#' @param iTERPSet S4 structure class ITERP
#' @param n.targ is the target columns in iTERPSet@Targets.spectra
#' @examples
#' plotEICs(object=iTERPSet, n.target = 1)


#t <- sapply(object@R,function(x) x[[target]][[1]][,ion])
#object@R[[sample]][[target]][[EICS]][,ion]

plotEICs <- function(object, n.target=1) {
  old.par <- par(no.readonly = T)
  n.targets <- ncol(object@Targets.spectra)
  par(mfrow=c(1, length(object@R[[1]][[n.targets]][[3]])))
  for(ion in 1:length(object@R[[1]][[n.target]][[3]])) {
    matplot(x=as.data.frame(lapply(object@R, function(x) as.numeric(rownames(x[[n.target]][["EICS"]])))),
            y=as.data.frame(lapply(object@R, function(x) x[[n.target]][["EICS"]][,ion])),type="l",
            xlab = "rt", ylab= "intensity",
            main = paste(names(object@R[[1]])[n.target],colnames(object@R[[1]][[n.target]]$EICS)[ion], sep="_")
    )
  }
  par(old.par)###preserves parameters
}

