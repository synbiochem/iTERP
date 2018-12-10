
#' getResults
#'
#' This function retrieves Area, Intensity and Extracted Ions Chromatograms for each Target
#' in real data. They are stores in slot R.
#'
#' @export getResults
#' @importFrom tcltk tk_choose.dir
#' @param n.ions number of ions we want to quantify for each targe
#' @examples
#' results <- getResults(object, type=c("I"))
#' results <- getResults(object, type=c("A"))



#object@R[[sample]][[target]]["EICS"]
getResults <- function(object, type=c("A","I")){
  r <- match.arg(type, c("A","I"), several.ok = FALSE)
  n.targets <- ncol(object@Targets.spectra)
  mydata <- NULL
  for (i in 1:n.targets){
    t1 <- lapply(object@R, "[[", i)
    t2 <- do.call(rbind,lapply(t1, "[[", r))
    colnames(t2) <- paste(colnames(object@Targets.spectra)[i], colnames(t2), sep="_")
    mydata[[i]] <- t2
  }
  results <- do.call(cbind,mydata)
  return(results)
}

