
#' Integrates TICs
#'
#' This function is aimed to establish the rt of target terpenoids in real data
#' @export integrateTIC
#' @param iTERPSet S4 structure class ITERP
#' @examples
#' ATIC <- integrateTIC(object=iTERPSet)
#'

integrateTIC <- function(object){
  A.TIC <- matrix(0, ncol=length(object@Target.comp.scan), nrow = length(object@TIC))
  for (Smp in 1:length(object@TIC)){
    x <- object@TIC[[Smp]]
    for (target in 1: length(object@Target.comp.scan)) {
      Kstart2 <- object@Target.comp.scan[target] - object@Targets.scan.window
      Kend2 <- object@Target.comp.scan[target] + object@Targets.scan.window
      A.TIC[Smp,target] <- sum(x[c(Kstart2:Kend2),"int"])
    }
  }
  rownames(A.TIC) <- names(object@TIC)
  colnames(A.TIC) <- colnames(object@Targets.spectra)
  return(A.TIC)
}
