.fit_clogit <- function (data, formula, more = NULL, summary_only = FALSE, ...) 
{
  message('maggies fit_clogit')
  if (!any(grepl(pattern = "^strata\\(.+\\)$", attr(terms(formula), 
                                                    "term.labels")))) {
    stop("No strata is provided, please make sure the formula includes a strata.")
  }
  m <- .clogit(formula, data = data, ...)
  if (summary_only) {
    m <- list(model = broom::tidy(m), 
              sl_ = attributes(data)$sl_, 
              ta_ = attributes(data)$ta_, more = more)
    class(m) <- c("fit_clogit_summary_only", "fit_clogit", 
                  class(m))
  }
  else {
    m <- list(model = m, sl_ = attributes(data)$sl_, ta_ = attributes(data)$ta_, 
              more = more)
    class(m) <- c("fit_clogit", class(m))
  }
  m
}





.clogit <- function (formula, data, weights, subset, na.action, 
                     method = c("exact", "approximate", "efron", "breslow"), ...) 
{
  message('maggies clogit')
  message(method)
  
  Call <- match.call()
  indx <- match(c("formula", "data"), names(Call), nomatch = 0)
  if (indx[1] == 0) 
    stop("A formula argument is required")
  mf <- Call[c(1, indx)]
  mf[[1L]] <- quote(stats::model.frame)
  mf$na.action <- "na.pass"
  nrows <- NROW(eval(mf, parent.frame()))
  coxcall <- Call
  coxcall[[1]] <- as.name(".coxph")
  newformula <- formula
  newformula[[2]] <- substitute(Surv(rep(1, nn), case), list(case = formula[[2]], 
                                                             nn = nrows))
  message('methods and formula')
  environment(newformula) <- environment(formula)
  coxcall$formula <- newformula
  method <- match.arg(method)
  coxcall$method <- switch(method, exact = "exact", efron = "efron", "breslow")
  if (method == "exact") {
    if (missing(data)) 
      temp <- terms(formula, specials = "cluster")
    else temp <- terms(formula, specials = "cluster", data = data)
    if (!is.null(attr(temp, "specials")$cluster) && method == 
        "exact") 
      stop("robust variance plus the exact method is not supported")
    if (!is.null(coxcall$weights)) {
      coxcall$weights <- NULL
      warning("weights ignored: not possible for the exact method")
    }
  }
  message('coxcall')
  coxcall <- eval(coxcall, sys.frame(sys.parent()))
  coxcall$userCall <- sys.call()
  class(coxcall) <- c("clogit", ".coxph")
  coxcall
}




.coxph <- function(formula, data, weights, subset, na.action,
                  init, control, ties= c("efron", "breslow", "exact"),
                  singular.ok =TRUE,  robust,
                  model=FALSE, x=FALSE, y=TRUE,  tt, method=ties, 
                  id, cluster, istate, statedata, nocenter=c(-1, 0, 1), ...) {
    message('maggies coxph')
    missing.ties <- missing(ties) & missing(method) #see later multistate sect
    ties <- match.arg(ties)
    Call <- match.call()
    ## We want to pass any ... args to coxph.control, but not pass things
    ##  like "dats=mydata" where someone just made a typo.  The use of ...
    ##  is simply to allow things like "eps=1e6" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs)) {
      controlargs <- names(formals(coxph.control)) #legal arg names
      indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
      if (any(indx==0L))
        stop(gettextf("Argument %s not matched", 
                      names(extraArgs)[indx==0L]), domain = NA)
    }
    if (missing(control)) control <- coxph.control(...) 
    
    # Move any cluster() term out of the formula, and make it an argument
    #  instead.  This makes everything easier.  But, I can only do that with
    #  a local copy, doing otherwise messes up future use of update() on
    #  the model object for a user stuck in "+ cluster()" mode.
    if (missing(formula)) stop("a formula argument is required")
    
    ss <- "cluster"
    if (is.list(formula))
      Terms <- if (missing(data)) terms(formula[[1]], specials=ss) else
        terms(formula[[1]], specials=ss, data=data)
    else Terms <- if (missing(data)) terms(formula, specials=ss) else
      terms(formula, specials=ss, data=data)
    
    tcl <- attr(Terms, 'specials')$cluster
    if (length(tcl) > 1) stop("a formula cannot have multiple cluster terms")
    
    if (length(tcl) > 0) { # there is one
      factors <- attr(Terms, 'factors')
      if (any(factors[tcl,] >1)) stop("cluster() cannot be in an interaction")
      if (attr(Terms, "response") ==0)
        stop("formula must have a Surv response")
      
      if (is.null(Call$cluster))
        Call$cluster <- attr(Terms, "variables")[[1+tcl]][[2]]
      else warning("cluster appears both in a formula and as an argument, formula term ignored")
      
      # [.terms is broken at least through R 4.1; use our
      #  local drop.special() function instead. 
      Terms <- drop.special(Terms, tcl)  
      formula <- Call$formula <- formula(Terms)
    }
    
    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    #  but don't evaluate it just yet
    indx <- match(c("formula", "data", "weights", "subset", "na.action",
                    "cluster", "id", "istate"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    tform <- Call[c(1,indx)]  # only keep the arguments we wanted
    tform[[1L]] <- quote(stats::model.frame)  # change the function called
    
    # if the formula is a list, do the first level of processing on it.
    if (is.list(formula)) {
      multiform <- TRUE
      dformula <- formula[[1]]   # the default formula for transitions   
      if (missing(statedata)) covlist <- parsecovar1(formula[-1])
      else {
        if (!inherits(statedata, "data.frame"))
          stop("statedata must be a data frame")
        if (is.null(statedata$state)) 
          stop("statedata data frame must contain a 'state' variable")
        covlist <- parsecovar1(formula[-1], names(statedata))
      }
      
      # create the master formula, used for model.frame
      # the term.labels + reformulate + environment trio is used in [.terms;
      #  if it's good enough for base R it's good enough for me
      tlab <- unlist(lapply(covlist$rhs, function(x) 
        attr(terms.formula(x$formula), "term.labels")))
      tlab <- c(attr(terms.formula(dformula), "term.labels"), tlab)
      newform <- reformulate(tlab, dformula[[2]])
      environment(newform) <- environment(dformula)
      formula <- newform
      tform$na.action <- na.pass  # defer any missing value work to later
    }
    else {
      multiform <- FALSE   # formula is not a list of expressions
      covlist <- NULL
      dformula <- formula
    }
    
    # add specials to the formula
    special <- c("strata", "tt", "frailty", "ridge", "pspline")
    tform$formula <- if(missing(data)) terms(formula, special) else
      terms(formula, special, data=data)
    
    # Make "tt" visible for coxph formulas, without making it visible elsewhere
    if (!is.null(attr(tform$formula, "specials")$tt)) {
      coxenv <- new.env(parent= environment(formula))
      assign("tt", function(x) x, envir=coxenv)
      environment(tform$formula) <- coxenv
    }
    
    # okay, now evaluate the formula
    mf <- eval(tform, parent.frame())
    Terms <- terms(mf)
    
    # Grab the response variable, and deal with Surv2 objects
    n <- nrow(mf)
    Y <- model.response(mf)
    isSurv2 <- inherits(Y, "Surv2")
    if (isSurv2) {
      # this is Surv2 style data
      # if there were any obs removed due to missing, remake the model frame
      if (length(attr(mf, "na.action"))) {
        tform$na.action <- na.pass
        mf <- eval.parent(tform)
      }
      if (!is.null(attr(Terms, "specials")$cluster))
        stop("cluster() cannot appear in the model statement")
      new <- surv2data(mf)
      mf <- new$mf
      istate <- new$istate
      id <- new$id
      Y <- new$y
      n <- nrow(mf)
    }       
    else {
      if (!is.Surv(Y)) stop("Response must be a survival object")
      id <- model.extract(mf, "id")
      istate <- model.extract(mf, "istate")
    }
    if (n==0) stop("No (non-missing) observations")
    if (length(id) >0) n.id <- length(unique(id))
    
    type <- attr(Y, "type")
    multi <- FALSE
    if (type=="mright" || type == "mcounting") multi <- TRUE
    else if (type!='right' && type!='counting')
      stop(paste("Cox model doesn't support \"", type,
                 "\" survival data", sep=''))
    data.n <- nrow(Y)   #remember this before any time transforms
    
    if (!multi && multiform)
      stop("formula is a list but the response is not multi-state")
    if (multi) {
      if (length(attr(Terms, "specials")$frailty) >0)
        stop("multi-state models do not currently support frailty terms")
      if (length(attr(Terms, "specials")$pspline) >0)
        stop("multi-state models do not currently support pspline terms")
      if (length(attr(Terms, "specials")$ridge) >0)
        stop("multi-state models do not currently support ridge penalties")
      if (!missing.ties) method <- ties <- "breslow"
    }
    
    if (control$timefix) Y <- aeqSurv(Y)
    if (length(attr(Terms, 'variables')) > 2) { # a ~1 formula has length 2
      ytemp <- innerterms(formula[1:2])
      suppressWarnings(z <- as.numeric(ytemp)) # are any of the elements numeric?
      ytemp <- ytemp[is.na(z)]  # toss numerics, e.g. Surv(t, 1-s)
      xtemp <- innerterms(formula[-2])
      if (any(!is.na(match(xtemp, ytemp))))
        warning("a variable appears on both the left and right sides of the formula")
    }
    
    # The time transform will expand the data frame mf.  To do this
    #  it needs Y and the strata.  Everything else (cluster, offset, weights)
    #  should be extracted after the transform
    #
    strats <- attr(Terms, "specials")$strata
    hasinteractions <- FALSE
    dropterms <- NULL
    if (length(strats)) {
      stemp <- untangle.specials(Terms, 'strata', 1)
      if (length(stemp$vars)==1) strata.keep <- mf[[stemp$vars]]
      else strata.keep <- strata(mf[,stemp$vars], shortlabel=TRUE)
      istrat <- as.integer(strata.keep)
      
      for (i in stemp$vars) {  #multiple strata terms are allowed
        # The factors attr has one row for each variable in the frame, one
        #   col for each term in the model.  Pick rows for each strata
        #   var, and find if it participates in any interactions.
        if (any(attr(Terms, 'order')[attr(Terms, "factors")[i,] >0] >1))
          hasinteractions <- TRUE  
      }
      if (!hasinteractions) dropterms <- stemp$terms 
    } else istrat <- NULL
    
    if (hasinteractions && multi)
      stop("multi-state coxph does not support strata*covariate interactions")
    
    
    timetrans <- attr(Terms, "specials")$tt
    if (missing(tt)) tt <- NULL
    if (length(timetrans)) {
      if (multi || isSurv2) stop("the tt() transform is not implemented for multi-state or Surv2 models")
      # begin tt() preprocessing
      timetrans <- untangle.specials(Terms, 'tt')
      ntrans <- length(timetrans$terms)
      
      if (is.null(tt)) {
        tt <- function(x, time, riskset, weights){ #default to O'Brien's logit rank
          obrien <- function(x) {
            r <- rank(x)
            (r-.5)/(.5+length(r)-r)
          }
          unlist(tapply(x, riskset, obrien))
        }
      }
      if (is.function(tt)) tt <- list(tt)  #single function becomes a list
      
      if (is.list(tt)) {
        if (any(!sapply(tt, is.function))) 
          stop("The tt argument must contain function or list of functions")
        if (length(tt) != ntrans) {
          if (length(tt) ==1) {
            temp <- vector("list", ntrans)
            for (i in 1:ntrans) temp[[i]] <- tt[[1]]
            tt <- temp
          }
          else stop("Wrong length for tt argument")
        }
      }
      else stop("The tt argument must contain a function or list of functions")
      
      if (ncol(Y)==2) {
        if (length(strats)==0) {
          sorted <- order(-Y[,1], Y[,2])
          newstrat <- rep.int(0L, nrow(Y))
          newstrat[1] <- 1L
        }
        else {
          sorted <- order(istrat, -Y[,1], Y[,2])
          #newstrat marks the first obs of each strata
          newstrat <-  as.integer(c(1, 1*(diff(istrat[sorted])!=0))) 
        }
        if (storage.mode(Y) != "double") storage.mode(Y) <- "double"
        counts <- .Call(Ccoxcount1, Y[sorted,], 
                        as.integer(newstrat))
        tindex <- sorted[counts$index]
      }
      else {
        if (length(strats)==0) {
          sort.end  <- order(-Y[,2], Y[,3])
          sort.start<- order(-Y[,1])
          newstrat  <- c(1L, rep(0, nrow(Y) -1))
        }
        else {
          sort.end  <- order(istrat, -Y[,2], Y[,3])
          sort.start<- order(istrat, -Y[,1])
          newstrat  <- c(1L, as.integer(diff(istrat[sort.end])!=0))
        }
        if (storage.mode(Y) != "double") storage.mode(Y) <- "double"
        counts <- .Call(Ccoxcount2, Y, 
                        as.integer(sort.start -1L),
                        as.integer(sort.end -1L), 
                        as.integer(newstrat))
        tindex <- counts$index
      }
      Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
      type <- 'right'  # new Y is right censored, even if the old was (start, stop]
      
      mf <- mf[tindex,]
      istrat <- rep(1:length(counts$nrisk), counts$nrisk)
      weights <- model.weights(mf)
      if (!is.null(weights) && any(!is.finite(weights)))
        stop("weights must be finite") 
      id <- model.extract(mf, "id")   # update the id and/or cluster, if present
      cluster <- model.extract(mf, "cluster")
      
      tcall <- attr(Terms, 'variables')[timetrans$terms+2]
      pvars <- attr(Terms, 'predvars')
      pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
      for (i in 1:ntrans) {
        newtt <- (tt[[i]])(mf[[timetrans$var[i]]], Y[,1], istrat, weights)
        mf[[timetrans$var[i]]] <- newtt
        nclass <- class(newtt)
        if (any(nclass %in% pmethod)) { # It has a makepredictcall method
          dummy <- as.call(list(as.name(class(newtt)[1]), tcall[[i]][[2]]))
          ptemp <- makepredictcall(newtt, dummy)
          pvars[[timetrans$terms[i]+2]] <- ptemp
        }
      }
      attr(Terms, "predvars") <- pvars
      # end tt() preprocessing
    }
    
    xlevels <- .getXlevels(Terms, mf)
    
    # grab the cluster, if present.  Using cluster() in a formula is no
    #  longer encouraged
    cluster <- model.extract(mf, "cluster")
    weights <- model.weights(mf)
    # The user can call with cluster, id, robust, or any combination
    # Default for robust: if cluster or any id with > 1 event or 
    #  any weights that are not 0 or 1, then TRUE
    # If only id, treat it as the cluster too
    has.cluster <- !(missing(cluster) || length(cluster)==0) 
    has.id <-      !(missing(id) || length(id)==0)
    has.rwt<-      (!is.null(weights) && any(weights != floor(weights)))
    #has.rwt<- FALSE  # we are rethinking this
    has.robust <-  (!missing(robust) && !is.null(robust))  # arg present
    if (has.id) id <- as.factor(id)
    
    if (missing(robust) || is.null(robust)) {
      if (has.cluster || has.rwt ||
          (has.id && (multi || anyDuplicated(id[Y[,ncol(Y)]==1]))))
        robust <- TRUE else robust <- FALSE
    }
    if (!is.logical(robust)) stop("robust must be TRUE/FALSE")
    
    if (has.cluster) {
      if (!robust) {
        warning("cluster specified with robust=FALSE, cluster ignored")
        ncluster <- 0
        clname <- NULL
      }
      else {
        if (is.factor(cluster)) {
          clname <- levels(cluster)
          cluster <- as.integer(cluster)
        } else {
          clname  <- sort(unique(cluster))
          cluster <- match(cluster, clname)
        }
        ncluster <- length(clname)
      }
    } else {
      if (robust && has.id) {
        # treat the id as both identifier and clustering
        clname <- levels(id)
        cluster <- as.integer(id)
        ncluster <- length(clname)
      }
      else {
        ncluster <- 0  # has neither
      }
    }
    
    # if the user said "robust", (time1,time2) data, and no cluster or
    #  id, complain about it
    if (robust && is.null(cluster)) {
      if (ncol(Y) ==2 || !has.robust) cluster <- seq.int(1, nrow(mf))
      else stop("one of cluster or id is needed") 
    }
    
    contrast.arg <- NULL  #due to shared code with model.matrix.coxph
    attr(Terms, "intercept") <- 1  # always have a baseline hazard
    
    if (multi) {
      # check for consistency of the states, and create a transition
      #  matrix
      if (length(id)==0) 
        stop("an id statement is required for multi-state models")
      
      mcheck <- survcheck2(Y, id, istate)
      # error messages here
      if (mcheck$flag["overlap"] > 0)
        stop("data set has overlapping intervals for one or more subjects")
      
      transitions <- mcheck$transitions
      istate <- mcheck$istate
      states <- mcheck$states
      
      #  build tmap, which has one row per term, one column per transition
      if (missing(statedata))
        covlist2 <- parsecovar2(covlist, NULL, dformula= dformula,
                                Terms, transitions, states)
      else covlist2 <- parsecovar2(covlist, statedata, dformula= dformula,
                                   Terms, transitions, states)
      tmap <- covlist2$tmap
      if (!is.null(covlist)) {
        # first vector will be true if there is at least 1 transition for which all
        #  covariates are present, second if there is at least 1 for which some are not
        good.tran <- bad.tran <- rep(FALSE, nrow(Y))  
        # We don't need to check interaction terms
        termname <- rownames(attr(Terms, 'factors'))
        trow <- (!is.na(match(rownames(tmap), termname)))
        
        # create a missing indicator for each term
        termiss <- matrix(0L, nrow(mf), ncol(mf))
        for (i in 1:ncol(mf)) {
          xx <- is.na(mf[[i]])
          if (is.matrix(xx)) termiss[,i] <- apply(xx, 1, any)
          else termiss[,i] <- xx
        }
        
        for (i in levels(istate)) {
          rindex <- which(istate ==i)
          j <- which(covlist2$mapid[,1] == match(i, states))  #possible transitions
          for (jcol in j) {
            k <- which(trow & tmap[,jcol] > 0)  # the terms involved in that 
            bad.tran[rindex] <- (bad.tran[rindex] | 
                                   apply(termiss[rindex, k, drop=FALSE], 1, any))
            good.tran[rindex] <- (good.tran[rindex] |
                                    apply(!termiss[rindex, k, drop=FALSE], 1, all))
          }
        }
        n.partially.used <- sum(good.tran & bad.tran & !is.na(Y))   
        omit <- (!good.tran & bad.tran) | is.na(Y)
        if (all(omit)) stop("all observations deleted due to missing values")
        temp <- setNames(seq(omit)[omit], attr(mf, "row.names")[omit])
        attr(temp, "class") <- "omit"
        mf <- mf[!omit,, drop=FALSE]
        attr(mf, "na.action") <- temp
        Y <- Y[!omit]
        id <- id[!omit]
        if (length(istate)) istate <- istate[!omit]  # istate can be NULL
      }
    }
    
    
    if (length(dropterms)) {
      Terms2 <- Terms[-dropterms]
      X <- model.matrix(Terms2, mf, constrasts.arg=contrast.arg)
      # we want to number the terms wrt the original model matrix
      temp <- attr(X, "assign")
      shift <- sort(dropterms)
      for (i in seq(along.with=shift))
        temp <- temp + 1*(shift[i] <= temp)
      attr(X, "assign") <- temp 
    }
    else X <- model.matrix(Terms, mf, contrasts.arg=contrast.arg)
    
    
    ## BOOKMARK
    print(summary(X[,'evi']))
    
    
    # drop the intercept after the fact, and also drop strata if necessary
    Xatt <- attributes(X) 
    if (hasinteractions) adrop <- c(0, untangle.specials(Terms, "strata")$terms)
    else adrop <- 0
    xdrop <- Xatt$assign %in% adrop  #columns to drop (always the intercept)
    X <- X[, !xdrop, drop=FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    attr(X, "contrasts") <- Xatt$contrasts
    offset <- model.offset(mf)
    if (is.null(offset) || all(offset==0)) {
      offset <- rep(0., nrow(mf))
      meanoffset <- 0
    } else if (any(!is.finite(offset) | !is.finite(exp(offset)))) 
      stop("offsets must lead to a finite risk score")
    else {
      meanoffset <- mean(offset)
      offset <- offset - meanoffset  # this can help stability of exp()
    }
    
    weights <- model.weights(mf)
    if (!is.null(weights) && any(!is.finite(weights)))
      stop("weights must be finite")   
    
    assign <- attrassign(X, Terms)
    contr.save <- attr(X, "contrasts")
    if (sum(Y[, ncol(Y)]) == 0) {
      # No events in the data!
      ncoef <- ncol(X)
      ctemp <- rep(NA, ncoef)
      names(ctemp) <- colnames(X)
      concordance= c(concordant=0, discordant=0, tied.x=0, tied.y=0, tied.xy=0,
                     concordance=NA, std=NA, timefix=FALSE)
      rval <- list(coefficients= ctemp,
                   var = matrix(0.0, ncoef, ncoef),
                   loglik=c(0,0),
                   score =0,
                   iter =0,
                   linear.predictors = offset,
                   residuals = rep(0.0, data.n),
                   means = colMeans(X), method=method,
                   n = data.n, nevent=0, terms=Terms, assign=assign,
                   concordance=concordance,  wald.test=0.0,
                   y = Y, call=Call)
      class(rval) <- "coxph"
      return(rval)
    }
    if (multi) {
      if (length(strats) >0) {
        # tmap starts with a "(Baseline)" row, which we want
        # strats is indexed off the data frame, which includes the response, so
        #  turns out to be correct for the remaining rows of tmap
        smap <- tmap[c(1L, strats),] 
        smap[-1,] <- ifelse(smap[-1,] >0, 1L, 0L)
        if (nrow(smap) > 2) {
          # multi state with more than 1 strata statement -- really unusual
          temp <- smap[-1,]
          if (!all(apply(temp, 2, function(x) all(x==0) || all(x==1)))) {
            # the hard case: some transitions use one strata variable, some
            #  transitions use another.  We need to keep them separate
            strata.keep <- mf[,strats]  # this will be a data frame
            istrat <- sapply(strata.keep, as.numeric)
          }
        }
      }
      else smap <- tmap[1,,drop=FALSE]
      cmap <- parsecovar3(tmap, colnames(X), attr(X, "assign"), covlist2$phbaseline)
      xstack <- stacker(cmap, smap, as.integer(istate), X, Y, strata=istrat,
                        states=states)
      
      rkeep <- unique(xstack$rindex)
      transitions <- survcheck2(Y[rkeep,], id[rkeep], istate[rkeep])$transitions
      
      Xsave <- X  # the originals may be needed later
      Ysave <- Y
      X <- xstack$X
      Y <- xstack$Y
      istrat <- xstack$strata
      if (length(offset)) offset <- offset[xstack$rindex]
      if (length(weights)) weights <- weights[xstack$rindex]
      if (length(cluster)) cluster <- cluster[xstack$rindex]
      t2 <- tmap[-c(1, strats),,drop=FALSE]   # remove the intercept row and strata rows
      r2 <- row(t2)[!duplicated(as.vector(t2)) & t2 !=0]
      c2 <- col(t2)[!duplicated(as.vector(t2)) & t2 !=0]
      a2 <- lapply(seq(along.with=r2), function(i) {cmap[assign[[r2[i]]], c2[i]]})
      # which elements are unique?  
      tab <- table(r2)
      count <- tab[r2]
      names(a2) <- ifelse(count==1, row.names(t2)[r2],
                          paste(row.names(t2)[r2], colnames(cmap)[c2], sep="_"))
      assign <- a2
    }
    
    # infinite covariates are not screened out by the na.omit routines
    #  But this needs to be done after the multi-X part
    if (!all(is.finite(X)))
      stop("data contains an infinite predictor")
    
    
    # init is checked after the final X matrix has been made
    if (missing(init)) init <- NULL
    else {
      if (length(init) != ncol(X)) stop("wrong length for init argument")
      temp <- X %*% init - sum(colMeans(X) * init) + offset
      # it's okay to have a few underflows, but if all of them are too
      #   small we get all zeros
      if (any(exp(temp) > .Machine$double.xmax) || all(exp(temp)==0))
        stop("initial values lead to overflow or underflow of the exp function")
    }
    
    pterms <- sapply(mf, inherits, 'coxph.penalty')
    if (any(pterms)) {
      pattr <- lapply(mf[pterms], attributes)
      pname <- names(pterms)[pterms]
      # 
      # Check the order of any penalty terms
      ord <- attr(Terms, "order")[match(pname, attr(Terms, 'term.labels'))]
      if (any(ord>1)) stop ('Penalty terms cannot be in an interaction')
      pcols <- assign[match(pname, names(assign))] 
      
      fit <- coxpenal.fit(X, Y, istrat, offset, init=init,
                          control,
                          weights=weights, method=method,
                          row.names(mf), pcols, pattr, assign, 
                          nocenter= nocenter)
    }
    else {
      rname <- row.names(mf)
      if (multi) rname <- rname[xstack$rindex]
      if( method=="breslow" || method =="efron") {
        if (grepl('right', type))  
          fit <- coxph.fit(X, Y, istrat, offset, init, control, 
                           weights=weights, method=method, 
                           rname, nocenter=nocenter)
        else  fit <- agreg.fit(X, Y, istrat, offset, init, control, 
                               weights=weights, method=method, 
                               rname, nocenter=nocenter)
      }
      else if (method=='exact') {
        if (type== "right")  
          fit <- coxexact.fit(X, Y, istrat, offset, init, control, 
                              weights=weights, method=method, 
                              rname, nocenter=nocenter)
        else fit <- agexact.fit(X, Y, istrat, offset, init, control, 
                                weights=weights, method=method, 
                                rname, nocenter=nocenter)
      }
      else stop(paste ("Unknown method to ties", method))
    }
    if (is.character(fit)) {
      fit <- list(fail=fit)
      class(fit) <- 'coxph'
    }
    else {
      if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
        vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
        msg <-paste("X matrix deemed to be singular; variable",
                    paste(vars, collapse=" "))
        if (!singular.ok) stop(msg)
        # else warning(msg)  # stop being chatty
      }
      fit$n <- data.n
      fit$nevent <- sum(Y[,ncol(Y)])
      if (length(id)>0) fit$n.id <- n.id
      fit$terms <- Terms
      fit$assign <- assign
      class(fit) <- fit$class
      fit$class <- NULL
      
      # don't compute a robust variance if there are no coefficients
      if (robust && !is.null(fit$coefficients) && !all(is.na(fit$coefficients))) {
        fit$naive.var <- fit$var
        # a little sneaky here: by calling resid before adding the
        #   na.action method, I avoid having missings re-inserted
        # I also make sure that it doesn't have to reconstruct X and Y
        fit2 <- c(fit, list(x=X, y=Y, weights=weights))
        if (length(istrat)) fit2$strata <- istrat
        if (length(cluster)) {
          temp <- residuals.coxph(fit2, type='dfbeta', collapse=cluster,
                                  weighted=TRUE)
          # get score for null model
          if (is.null(init))
            fit2$linear.predictors <- 0*fit$linear.predictors
          else fit2$linear.predictors <- c(X %*% init)
          temp0 <- residuals.coxph(fit2, type='score', collapse=cluster,
                                   weighted=TRUE)
        }
        else {
          temp <- residuals.coxph(fit2, type='dfbeta', weighted=TRUE)
          fit2$linear.predictors <- 0*fit$linear.predictors
          temp0 <- residuals.coxph(fit2, type='score', weighted=TRUE)
        }
        fit$var <- t(temp) %*% temp
        u <- apply(as.matrix(temp0), 2, sum)
        fit$rscore <- coxph.wtest(t(temp0)%*%temp0, u, control$toler.chol)$test
      }
      
      #Wald test
      if (length(fit$coefficients) && is.null(fit$wald.test)) {  
        #not for intercept only models, or if test is already done
        nabeta <- !is.na(fit$coefficients)
        # The init vector might be longer than the betas, for a sparse term
        if (is.null(init)) temp <- fit$coefficients[nabeta]
        else temp <- (fit$coefficients - 
                        init[1:length(fit$coefficients)])[nabeta]
        fit$wald.test <-  coxph.wtest(fit$var[nabeta,nabeta], temp,
                                      control$toler.chol)$test
      }
      
      # Concordance.  Done here so that we can use cluster if it is present
      # The returned value is a subset of the full result, partly because it
      #  is all we need, but more for backward compatability with survConcordance.fit
      if (length(cluster))
        temp <- concordancefit(Y, fit$linear.predictors, istrat, weights,
                               cluster=cluster, reverse=TRUE,
                               timefix= FALSE)
      else temp <- concordancefit(Y, fit$linear.predictors, istrat, weights,
                                  reverse=TRUE, timefix= FALSE)
      if (is.matrix(temp$count))
        fit$concordance <- c(colSums(temp$count), concordance=temp$concordance,
                             std=sqrt(temp$var))
      else fit$concordance <- c(temp$count, concordance=temp$concordance, 
                                std=sqrt(temp$var))
      
      na.action <- attr(mf, "na.action")
      if (length(na.action)) fit$na.action <- na.action
      if (model) {
        if (length(timetrans)) {
          stop("'model=TRUE' not supported for models with tt terms")
        }
        fit$model <- mf
      }
      if (x)  {
        if (multi) fit$x <- Xsave else fit$x <- X
        if (length(timetrans)) fit$strata <- istrat
        else if (length(strats)) fit$strata <- strata.keep
      }
      if (y)  {
        if (multi) fit$y <- Ysave else fit$y <- Y
      }
      fit$timefix <- control$timefix  # remember this option
    }
    if (!is.null(weights) && any(weights!=1)) fit$weights <- weights
    if (multi) {
      fit$transitions <- transitions
      fit$states <- states
      fit$cmap <- cmap
      fit$smap <- smap   # why not 'stratamap'?  Confusion with fit$strata
      nonzero <- which(colSums(cmap)!=0)
      fit$rmap <- cbind(row=xstack$rindex, transition= nonzero[xstack$transition])
      
      # add a suffix to each coefficent name.  Those that map to multiple transitions
      #  get the first transition they map to
      single <- apply(cmap, 1, function(x) all(x %in% c(0, max(x)))) #only 1 coef
      cindx <- col(cmap)[match(1:length(fit$coefficients), cmap)]
      rindx <- row(cmap)[match(1:length(fit$coefficients), cmap)]
      suffix <- ifelse(single[rindx], "", paste0("_", colnames(cmap)[cindx]))
      newname <- paste0(names(fit$coefficients), suffix)
      if (any(covlist2$phbaseline > 0)) {
        # for proportional baselines, use a better name
        base  <- colnames(tmap)[covlist2$phbaseline]
        child <- colnames(tmap)[which(covlist2$phbaseline >0)]
        indx <- 1 + length(newname) - length(base):1 # coefs are the last ones
        newname[indx] <-  paste0("ph(", child, "/", base, ")")
        phrow <- apply(cmap, 1, function(x) all(x[x>0] %in% indx))
        matcoef <- cmap[!phrow,,drop=FALSE ] # ph() terms exluded 
      }
      else matcoef <- cmap   
      names(fit$coefficients) <- newname
      
      if (FALSE) { 
        # an idea that was tried, then paused: make the linear predictors
        # and residuals into matrices with one column per transition
        matcoef[matcoef>0] <- fit$coefficients[matcoef]
        temp <- Xsave %*% matcoef
        colnames(temp) <- colnames(cmap)
        fit$linear.predictors <- temp
        
        temp <- matrix(0., nrow=nrow(Xsave), ncol=ncol(fit$cmap))
        temp[cbind(xstack$rindex, xstack$transition)] <- fit$residuals
        # if there are any transitions with no covariates, residuals have not
        #  yet been calculated for those.
        if (any(colSums(cmap) ==0)) {
          from.state <- as.numeric(sub(":.*$", "", colnames(cmap)))
          to.state   <- as.numeric(sub("^.*:", "", colnames(cmap)))
          # warning("no covariate residuals not filled in")
        }
        fit$residuals <- temp
      }
      class(fit) <- c("coxphms", class(fit))
    }
    names(fit$means) <- names(fit$coefficients)
    
    fit$formula <- formula(Terms)
    if (length(xlevels) >0) fit$xlevels <- xlevels
    fit$contrasts <- contr.save
    if (meanoffset !=0) fit$linear.predictors <- fit$linear.predictors + meanoffset
    if (x & any(offset !=0)) fit$offset <- offset
    
    fit$call <- Call
    fit
}















vcov.coxph <- function (object, complete=TRUE, ...) {
  # conform to the standard vcov results
  vmat <- object$var
  vname <- names(object$coefficients)
  dimnames(vmat) <- list(vname, vname)
  if (!complete && any(is.na(coef(object)))) {
    keep <- !is.na(coef(object))
    vmat[keep, keep, drop=FALSE]
  }
  else vmat
}

vcov.survreg<-function (object, complete=TRUE, ...) {
  if (!complete && any(is.na(coef(object)))) {
    keep <- !is.na(coef(object))
    vv <- object$var[keep, keep, drop=FALSE]
    vname <- names(coef(object))[keep]
  }
  else {
    vv <- object$var
    vname <- names(coef(object))   # add dimnames
  }
  extra <- ncol(vv) - length(vname)
  if (extra ==1) vname <- c(vname, "Log(scale)")
  else if(extra >1) 
    vname <- c(vname, paste("Log(scale[", names(object$scale), "])", sep=''))
  dimnames(vv) <- list(vname, vname)
  vv
}

# The extractAIC methods for coxph and survreg objects are defined
#  in the stats package.  Don't reprise them here.
extractAIC.coxph.penal<- function(fit,scale,k=2,...){
  edf<-sum(fit$df)
  loglik <- fit$loglik[length(fit$loglik)]
  c(edf, -2 * loglik + k * edf)
}

extractAIC.coxph.null <- function(fit, scale, k=2, ...) {
  c(0, -2*fit$loglik[1])
}

labels.survreg <- function(object, ...) attr(object$terms, "term.labels")


# This function is just like all.vars -- except that it does not recur
#  on the $ sign, it follows both arguments of +, *, - and : in order to
#  track formulas, all arguments of Surv, and only the first of things 
#  like ns().  And - it works only on formulas.
# This is used to generate a warning in coxph if the same variable is used
#  on both sides, so perfection is not required of the function.
# Changed from terms.inner to innerterms at CRAN's request, as the former
#  created a false positive as an undocumented method for the terms generic
innerterms <- function(x) {
  if (inherits(x, "formula")) {
    if (length(x) ==3) c(innerterms(x[[2]]), innerterms(x[[3]]))
    else innerterms(x[[2]])
  }
  else if (inherits(x, "call") && 
           (x[[1]] != as.name("$") && x[[1]] != as.name("["))) {
    if (x[[1]] == '+' || x[[1]]== '*' || x[[1]] == '-' || x[[1]] ==':') {
      # terms in a model equation, unary minus only has one argument
      if (length(x)==3) c(innerterms(x[[2]]), innerterms(x[[3]]))
      else innerterms(x[[2]])
    }
    else if (x[[1]] == as.name("Surv"))
      unlist(lapply(x[-1], innerterms))
    else if (length(x) ==2) innerterms(x[[2]])
    else character(0)
  }
  else(deparse(x))
}

# If a subject had (start, stop) observations of (1,2) (2,10) (10,15) (20,25),
#  say, code often wants to distiguish intervals that are "real" censoring
#  from a simple split due to a time dependent covariate.
# To support users who like the 'extended Kaplan-Meier' , i.e. a person can
#  switch curves midstream, we need to define someone as censored from the
#  first curve and an entry to the second whenever such a switch occurs.
# This routine returns 1*(first in a sequence) + 2*(last in a sequence),
#  which for the above is 1,0,2,3.  This assumes no overlapping intervals.

survflag <- function(y, id, group) {
  if (!inherits(y, "Surv")) stop("y must be a Surv object")
  if (nrow(y) != length(id)) stop("length mismatch")
  if (ncol(y) != 3) stop("y needs to be of (tstart, tstop) form")
  
  n <- nrow(y)
  if (missing(group))
    indx <- order(id, y[,2])  # sort the data by time within id
  else indx <- order(group, id, y[,2])
  y2 <- y[indx,]
  id2 <- id[indx]
  
  if (missing(group)) newid <- (id2[-n] != id2[-1])
  else {
    group2 <-as.numeric(group)[indx]  #normally group is a factor
    newid <- ((id2[-n] != id2[-1]) | (group2[-n] != group2[-1]))
  }       
  gap <-  (y2[-n,2] < y2[-1,1]) 
  
  flag <- unname(1L*c(TRUE, newid | gap) + 2L*c(newid | gap, TRUE))
  flag[indx] <- flag   # return it to data order
  flag
}

# Dummy methods, to create an informative error message
coef.survfit <- function(object, ...) 
  stop("coef method not applicable for survfit objects")
vcov.survfit <- function(object, ...) 
  stop("vcov method not applicable for survfit objects")
confint.survfit <- function(object, ...)
  stop(paste("confint method not defined for survfit objects," ,
             "use quantile for confidence intervals of the median survival"))













# This routine fits right censored data when the method is 
#  "exact".  The most common use for this option is matched
#   case-control data.
coxexact.fit <- function(x, y, strata, offset, init, control,
                         weights, method, rownames,
                         resid=TRUE, nocenter=NULL)
{
  if (!is.matrix(x)) stop("Invalid formula for cox fitting function")
  if (!is.null(weights) && any(weights!=1))
    stop("Case weights are not supported for the exact method")
  n <- nrow(x)
  nvar <- ncol(x)
  
  # The risk set addition in the C-code, which is the critically slow
  # part of the calculations, expects to have the data in sorted order:
  #   (large to small times) within strata
  if (length(strata)==0) {
    sorted <- order(-y[,1])
    newstrat <- as.integer(rep(0,n))
  }
  else {
    sorted <- order(strata, -y[,1])
    strata <- (as.numeric(strata))[sorted]
    newstrat <- as.integer(c(1, 1*(diff(strata)!=0)))
  }
  y <- y[sorted,]
  if (is.null(offset)) offset <- rep(0.,n)
  else offset <- offset[sorted]
  
  if (nvar==0) {
    # A special case: Null model.  Trick the C code, which requires
    #   at least one variable, by creating one and then doing 0
    #   iterations at beta=0
    x <- matrix(1:n, ncol=1)
    init <- NULL
    maxiter <- 0
    nullmodel <- TRUE
    nvar <- 1
  }
  else {
    maxiter <- control$iter.max
    nullmodel <- FALSE
  }
  
  if (!is.null(init)) {
    if (length(init) != nvar) stop("Wrong length for inital values")
  }
  else init <- rep(0.,nvar)
  
  # Prescale the data set to improve numerical accuracy.
  #  We will undo the scaling before finishing up.
  newx <- scale(x[sorted,])
  rescale <- attr(newx, "scaled:scale")
  means   <- attr(newx, "scaled:center")
  if (!is.null(nocenter)){
    zero.one <- apply(x, 2, function(z) all(z %in% nocenter))
    for (i in which(zero.one)) {
      newx[,i] <- x[sorted,i]
      rescale[i] <- 1.0
      means[i]   <- 0.0
    }
  }
  cfit <- .Call(Ccoxexact, 
                as.integer(maxiter),
                as.double(y),  # integer data?  Just in case.
                newx,
                as.double(offset),
                as.integer(newstrat),
                as.double(init*rescale),
                as.double(control$eps),
                as.double(control$toler.chol)
  )
  if (nullmodel) {
    score <- exp(offset[sorted])
    cxres <- .C(Ccoxmart2,
                as.integer(n),
                as.double(y[,1]),
                as.integer(y[,2]),
                newstrat,
                score,
                rep(1.0, n),  #weights
                resid=double(n))
    resid <- double(n)
    resid[sorted] <- cxres$resid
    names(resid) <- rownames
    return( list(loglik = cfit$loglik[1],
                 linear.predictors = offset,
                 residuals = resid,
                 method = method,
                 class = c("coxph.null", "coxph")))
  }
  
  loglik <- cfit$loglik[1:2]  #these are packed into one vector
  sctest <- cfit$loglik[3]
  iter <- cfit$loglik[5]
  flag <- cfit$loglik[4]
  var <- matrix(cfit$imat,nvar,nvar)
  coef <- cfit$coef
  
  if (flag < nvar) which.sing <- diag(var)==0
  else which.sing <- rep(FALSE,nvar)
  
  infs <- abs(cfit$u %*% var)
  if (control$iter.max >1) {
    if (flag == 1000)
      warning("Ran out of iterations and did not converge")
    else {
      infs <- ((infs > control$eps) & 
                 infs > control$toler.inf*abs(coef))
      if (any(infs))
        warning(paste("Loglik converged before variable ",
                      paste((1:nvar)[infs],collapse=","),
                      "; beta may be infinite. "))
    }
  }
  
  names(coef) <- dimnames(x)[[2]]
  lp  <- newx %*% coef + offset 
  score <- as.double(exp(lp))
  
  if (resid) {
    # Compute the residuals
    cxres <- .C(Ccoxmart2,
                as.integer(n),
                as.double(y[,1]),
                as.integer(y[,2]),
                newstrat,
                score,
                rep(1.0, n),  #weights
                resid=double(n))
    resid <- double(n)
    resid[sorted] <- cxres$resid
    names(resid) <- rownames
    coef[which.sing] <- NA
    lp.unsort <- double(n)
    lp.unsort[sorted] <- lp
    
    scmat <- diag(1/rescale, nvar,nvar)
    rval <-  list(coefficients  = coef/rescale,
                  var    = scmat %*% var %*% scmat,
                  loglik = loglik,
                  score  = sctest,
                  iter   = iter,
                  linear.predictors = lp.unsort,
                  residuals = resid,
                  means = means,
                  method = method,
                  class = 'coxph')
  } else {
    rval <-  list(coefficients  = coef/rescale,
                  var    = scmat %*% var %*% scmat,
                  loglik = loglik,
                  score  = sctest,
                  iter   = iter,
                  linear.predictors = lp.unsort,
                  means = means,
                  method = method,
                  class = 'coxph')
  }
  rval
}