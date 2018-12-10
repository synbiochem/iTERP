#'  NIST2Targets
#'
#' Reads .MSP files for a directory and saves in to a target matrix
#'
#' aaaaaaaaaaaaaaaaaaaaaaaaaaaa
#'
#' @export NIST2targets
#' @export varEntryDialog
#' @exportClass iTERP
#' @import tcltk
#' @importFrom OrgMassSpecR ReadMspFile
#' @usage iTERPSet <- NIST2targets()



iTERPSet <- setClass("iTERP", slots=c(Targets.spectra="data.frame",
                                      Targets.rt.input="list",
                                      Targets.ion.input="list",
                                      Sample.to.locate.targets="character",
                                      Targets.scan.window = "numeric",
                                      Targets.rt ="numeric",
                                      Targets.rt.vectors ="list",
                                      Target.comp.scan = "numeric",
                                      Sample.name = "character",
                                      EICs= "list",
                                      As = "list",
                                      Is = "list",
                                      TICs= "list",
                                      TIC = "list",
                                      R = "list",
                                      IS.curve="list"))


varEntryDialog <- function(vars,
                           labels = vars,
                           fun = rep(list(as.character), length(vars)),
                           title = 'Variable Entry',
                           prompt = NULL) 
{
  
  stopifnot(length(vars) == length(labels), length(labels) == length(fun))
  
  # Create a variable to keep track of the state of the dialog window:
  # done = 0; If the window is active
  # done = 1; If the window has been closed using the OK button
  # done = 2; If the window has been closed using the Cancel button or destroyed
  done <- tclVar(0)
  
  tt <- tktoplevel()
  tkwm.title(tt, title)
  entries <- list()
  tclvars <- list()
  
  # Capture the event "Destroy" (e.g. Alt-F4 in Windows) and when this happens,
  # assign 2 to done.
  tkbind(tt,"<Destroy>",function() tclvalue(done)<-2)
  
  for(i in seq_along(vars)) {
    tclvars[[i]] <- tclVar("")
    entries[[i]] <- tkentry(tt, textvariable=tclvars[[i]])
  }
  
  doneVal <- as.integer(tclvalue(done))
  results <- list()
  
  reset <- function() {
    for(i in seq_along(entries)) {
      tclvalue(tclvars[[i]]) <<- ""
    }
  }
  reset.but <- tkbutton(tt, text="Reset", command=reset)
  
  cancel <- function() {
    tclvalue(done) <- 2
  }
  cancel.but <- tkbutton(tt, text='Cancel', command=cancel)
  
  submit <- function() {
    for(i in seq_along(vars)) {
      tryCatch( {
        results[[vars[[i]]]] <<- fun[[i]](tclvalue(tclvars[[i]]))
        tclvalue(done) <- 1
      },
      error = function(e) { tkmessageBox(message=geterrmessage()) },
      finally = { }
      )
    }
  }
  submit.but <- tkbutton(tt, text="Submit", command=submit)
  
  if(!is.null(prompt)) {
    tkgrid(tklabel(tt,text=prompt), columnspan=3, pady=10)
  }
  
  for(i in seq_along(vars)) {
    tkgrid(tklabel(tt, text=labels[i]), entries[[i]], pady=10, padx=10, columnspan=4)
  }
  
  tkgrid(submit.but, cancel.but, reset.but, pady=10, padx=10, columnspan=3)
  tkfocus(tt)
  
  # Do not proceed with the following code until the variable done is non-zero.
  #   (But other processes can still run, i.e. the system is not frozen.)
  tkwait.variable(done)
  
  if(tclvalue(done) != 1) {
    results <- NULL
  }
  
  tkdestroy(tt)
  return(results)
}

rel.int <- function(x){(x/(max(x)))*100}

readNIST <- function(file.Name) {
    conn <- file(file.Name,open="r")
    linn <-readLines(conn)
    n.peaks <- as.numeric(unlist(strsplit(linn[grep("Num Peaks:",linn)],split = ":"))[2])
    MS <- OrgMassSpecR::ReadMspFile(file=file.Name,skip = grep("Num Peaks:",linn))
    MS$intensity <- rel.int(MS$intensity)
    close(conn)
    return(MS)
}

NIST2targets <- function(rt.sel=F, ion.sel=F) {
  nom <- tk_choose.files(caption = "Select .MSP files to targets",
                         multi = TRUE)

  n.target <- length(nom)
  spTarg <- matrix(0, ncol=n.target, nrow=600)
  col.name <- NULL
  for(i in 1:n.target){
    file.Name <- nom[i]
    columna <- unlist(strsplit(file.Name, split = "/"))
    col.name[i] <- gsub(".MSP","",columna[length(columna)])
    MS <- readNIST(file.Name)
    spTarg[MS$mz,i] <- MS$intensity
  }
  spTarg <- data.frame(spTarg)
  colnames(spTarg) <- col.name
  rownames(spTarg) <- c(1:nrow(spTarg))

  object <- new("iTERP", Targets.spectra = spTarg)
  targets <- colnames(spTarg)
  
  if(rt.sel==T){
    rt <- varEntryDialog(vars=targets, labels = targets,
                         title = "Enter rt")
    rt <- lapply(rt, function(x) as.numeric(x))
  } else{
    rt <- as.list(rep(NA, times=length(targets)))
  }
  
  if(ion.sel==T){
    ion <- varEntryDialog(vars=targets, labels = targets,
                         title = "Enter selective ion")
    ion <- lapply(ion, function(x) as.numeric(x))
  } else{
    ion <- as.list(rep(NA, times=length(targets)))
  }

  object@Targets.rt.input <- rt
  object@Targets.ion.input <- ion
  
  return(object)
}

