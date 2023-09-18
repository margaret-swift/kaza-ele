# noodle.R
# Created XX DATE XX
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/noodle.R')
source('02_scripts/utilities.R')
# file = list.files(here::here('01_data', 'foldername', 'raw'), full.names = TRUE)


# ******************************************************************************
#                             MESSING AROUND WITH OOP
# ******************************************************************************
# https://adv-r.hadley.nz/s4.html

# ******************************************************************************
# __init__

### ======== ELEPHANT CLASS ======== ###

#' DEFINITION for the Elephant Class
setClass("Elephant",
  slots=list(ID=   "numeric",
             SEX=  "character",
             X=    "numeric",
             Y=    "numeric",
             HISTORY= "list",
             PAR1= "numeric"))

#' HELPER for the Elephant class
#' - takes in the SEX slot and transforms it into a uniform M or F,
#' - pulls hyperparameters based on the SEX, and uses these to randomly
#'     generate various parameters that we will use in the RWs. HPs should
#'     come from the literature.
ele <- function(ID, SEX = "F") {
  
  # set hyper-parameters based on sex
  hypers <-  list("M" = list(mu=0, sd=1),
                  "F" = list(mu=5, sd=2.5))
  SEX <- toupper(substr(SEX, 0, 1))
  if (SEX %in% c("M", "F")) { 
    pars <- hypers[[SEX]]
    hmu = pars$mu
    hsd = pars$sd
  } else { warning('SEX should be M or F') }
  
  # generate parameters based on hyper-parameters
  PAR1 = rnorm(1, hmu, hsd)
  
  # randomly generate initial coordinates
  X = rnorm(1, 0, 1)
  Y = rnorm(1, 0, 1)
  
  new("Elephant", 
      ID = ID, 
      SEX = SEX,
      X=X,
      Y=Y,
      HISTORY = list(list(X=X, Y=Y),
                     list(X=X, Y=Y)),
      PAR1=PAR1)
}
e=ele(123, 'f')
e


#' VIEWER for Elephant class
#' Defining a function to display object details
setMethod("show", "Elephant",
  function(object){
    
    makeLoc <- function(x, y) paste0(round(x, 2), ', ', round(y, 2))
    cur = makeLoc(object@X, object@Y)
    hist = lapply(object@HISTORY, function(e) makeLoc(e$X, e$Y))[1:3]
    hist = paste(paste(hist, collapse='\n     '), '\n    ...\n')
    cat(is(object)[[1]], '\n',
        "  ID: ", object@ID, "\n",
        "  SEX: ", object@SEX, "\n",
        "  PAR1: ", object@PAR1, "\n",
        "  LOCATION: ", cur, '\n',
        "  HISTORY:\n    ", hist, '\n')
  }
)

#' ACCESSORS to Elephant parameters
par <- function(x) x@PAR1
xloc <- function(x) x@X
yloc <- function(x) x@Y
sex <- function(x) x@SEX
id <- function(x) x@ID
history <- function(x) x@HISTORY

#' SETTERS to Elephant parameters
setGeneric("par<-", function(x, value) standardGeneric("par<-"))
setMethod("par<-", "Elephant", function(x, value) {obj@PAR1 <- value; obj})
setGeneric("x<-", function(obj, value) standardGeneric("x<-"))
setMethod("x<-", "Elephant", function(obj, value) {obj@Y <- value; obj})
setGeneric("y<-", function(obj, value) standardGeneric("y<-"))
setMethod("y<-", "Elephant", function(obj, value) {obj@Y <- value; obj})


### ================================ ###


### ======== OTHER METHODS ======== ###


#' Generic to take a step
setGeneric('takeStep', function(x) standardGeneric("takeStep"))
setMethod("takeStep", "Elephant", 
  function(x) {
    print(x@PAR1) 
  }
)


# ******************************************************************************
# __main__

# TESTING

set.seed(2424)
nrep=5e3
eles <- data.frame(ID = 1:nrep, SEX = "", PAR1 = 0)
for (i in 1:nrep) {
  si <- ifelse(round(runif(1)), "F", "M")
  ei <- ele(ID = i, SEX = si)
  eles$SEX[i] <- si
  eles$PAR1[i] <- par(ei)
}

ggplot() + geom_histogram(data=eles, 
                          aes(fill=SEX, x=PAR1), 
                          alpha=0.7, position="identity", bins=100)




