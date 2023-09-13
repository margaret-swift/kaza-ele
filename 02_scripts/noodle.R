# noodle.R
# Created XX DATE XX
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

# here::i_am('02_scripts/noodle.R')
# source('02_scripts/utilities.R')
# file = list.files(here::here('01_data', 'foldername', 'raw'), 
#                  full.names = TRUE)

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************



# ******************************************************************************
#                             MESSING AROUND WITH OOP
# ******************************************************************************


# ******************************************************************************
# __init__
ele <- setClass("Elephant", 
                slots=list(ID="numeric",
                           SEX="character",
                           X='numeric',
                           Y='numeric',
                           PAR1="numeric"))
hypers <-  list("M" = list(mu=0, sd=1),
                "F" = list(mu=5, sd=2.5))


### ======== INITIALIZING CLASSES ======== ###
setMethod("initialize", "Elephant", function(.Object, ...) {
  .Object <- callNextMethod()
  .Object@SEX <- toupper(substr(.Object@SEX, 0, 1))
  
  # set hyper-parameters based on sex
  if (.Object@SEX %in% c("M", "F")) { 
    pars <- hypers[[.Object@SEX]]
    hmu = pars$mu
    hsd = pars$sd
  } else {
    warning('SEX should be M or F')
  }
  
  # generate parameters based on hyper-parameters
  .Object@PAR1 = rnorm(1, hmu, hsd)
  
  # randomly generate initial coordinates
  .Object@X = rnorm(1, 0, 1)
  .Object@Y = rnorm(1, 0, 1)
  .Object
})





### ======== OTHER METHODS ======== ###

# Defining a function to display object details
setMethod("show", "Elephant",
  function(object){
    cat("ID:  ", object@ID, "\n")
    cat("SEX: ", object@SEX, "\n")
    cat("PAR1: ", object@PAR1, "\n")
    cat('LOCATION: ', paste0(object@X, ', ', object@Y), '\n')
  }
)
# updates
setMethod('update', "Elephant",
  function(object){
    object@PAR1 <- object@PAR1 + 1
    object
  }
)

# setting generic function calls
setGeneric('par', function(x) standardGeneric("par"))
setGeneric("par<-", function(x, value) standardGeneric("par<-"))
setMethod("par", "Elephant", function(x) x@PAR1)
setMethod("par<-", "Elephant", function(x, value) {x@PAR1 <- value; x})

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









