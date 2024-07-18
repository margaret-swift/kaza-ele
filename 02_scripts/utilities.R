# utilities.R
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                    LIBRARY LOADING & DEFINE GLOBAL OBJECTS
# ******************************************************************************
message("\nLoading all utility functions and parameters from utilities.R\n")

message('Loading base packages for all projects...')
pacman::p_load(tidyverse, patchwork, sf, lubridate, tictoc, here)
sf_use_s2(FALSE) #fix errors for spherical geometry
message('   ...Base packages loaded.')

message("Setting base objects...")
datadir <- here::here("01_data")
scriptdir <- here::here("02_scripts")
outdir <- here::here("03_output")
crs_LL <- 'EPSG:4326'
crs_UTM<- 'EPSG:32734'
message("   ...Base objects loaded.")

# Plot Params
message("Setting plot params...")
plot.theme <- theme( text = element_text(size=18), 
                     title = element_text(size=20),
                     axis.text.x = element_text(angle=45, vjust=0.5))
big.theme = theme(text=element_text(size=16), title = element_text(size=20))
blank.theme <- theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
transp.theme <- theme(
  panel.background = element_rect(fill='transparent'),
  plot.background = element_rect(fill='transparent', color=NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.background = element_rect(fill='transparent'),
  legend.box.background = element_rect(fill='transparent')
)
# scalebar <- ggspatial::annotation_scale( pad_x = unit(0.05, "in"), 
#                                          pad_y = unit(0.05, "in"),
#                                          text_face="bold") 
message("   ...Plot params loaded.")

# ******************************************************************************
#                             BASE UTILITY FUNCTIONS
# ******************************************************************************
message("Loading functions...")

as.c <- as.character
as.n <- as.numeric
as.v <- as.vector
nog  <- st_drop_geometry
len  <- length

removePunct <- function(v) gsub('\\)|\\(|\\/|\\.|\\?', '', v)
replaceSpace <- function(v, r="_") gsub(' ', r, v)
fixNames <- function(d) {
  names(d) <- replaceSpace(removePunct(toupper(names(d))))
  d
}

# functions to set user-friendly datapath options
metapath <- function(...) .metapath(...)
rawpath  <- function(...) .rawpath(...)
procpath <- function(...) .procpath(...)
outpath  <- function(...) .outpath(...)
setOutPath <- makeOutPath <- function(name, 
                                      create=FALSE, 
                                      verbose=TRUE) {
  if (length(name)>1) name <- paste(name, collapse="/")
  path <- here::here(outdir, name)
  if (!dir.exists(path)) {
    if (create) {
      message('\nOutput directory "', name, '" created. ')
      dir.create(path , showWarnings = FALSE)
    } else {
      message('\nOutput directory "', name, '" not found.\n',
        'To create new directory, use "create=TRUE".\n',
        'Or, check spelling against existing paths: ')
      cat(paste(' -', list.files(outdir), collapse="\n"))
    }
  } 
  if (dir.exists(path)){
    cat('rebasing outpath() paths to ./03_output/', name, '\n',
        'Existing files in directory:\n', sep='')
    .outpath <<- function(...) here::here(path, ...)
    if (verbose) { message(paste(' -', list.files(path), collapse="\n")) }
  }
}
setDataPaths <- function(dataname, verbose=TRUE, print.length=10) {
  dp <- here::here(datadir, dataname)
  if (!dir.exists(dp)) {
    warning('\n  Data directory "', dataname,
            '" not found. Did you mean one of these in /0_data/? ',
            paste("\n     -", list.files(datadir)))
  } else {
    .metapath <<- function(...) here::here(dp, 'meta', ...)
    .rawpath  <<- function(...) here::here(dp, 'raw', ...)
    .procpath <<- function(...) here::here(dp, 'processed', ...)

    if (verbose) {
      # function to easily print path names and files within
      printPathFiles <- function(slug) {
        # print directory name
        shortslug <- substr(slug, 0, 4)
        fname <- file.path(gsub(".*01_data", "", dp), slug)
        fprint<- ifelse(str_count(shortslug) == 4,'path() <-','path()  <-' )
        cat('  ', paste0(shortslug, fprint), fname)

        # print files
        files <- list.files(file.path(dp, slug))
        L=min(length(files), print.length)
        message(paste("\n     -", files[1:L]))
        if (length(files) > print.length) message('     ...')
      }
      cat('rebasing data grabber function paths to...\n')
      printPathFiles('meta')
      printPathFiles('raw')
      printPathFiles('processed')
    }
  }
}
quickload <- function(dnames = c('linear_features', 'boundaries')) {
  loadme <- function(dn) {
    datalist <- list.files(datadir)
    if ( dn %in% datalist) {
      fname <- paste0(dn, ".rdata")
      datapath <- here::here(datadir, dn, "processed")
      path <- here::here(datapath, fname)
      if (file.exists(path)) { 
        load(file=path, envir=globalenv())
        message('Loaded ', toupper(dn), ' data into global environment.')
      } else { 
        message("WARNING: ", fname, " does not exist. Available files: ")
        print(datapath)
        cat(paste(' -', list.files(datapath), collapse="\n" ))
      }
    } else {
      message("WARNING: ", toupper(dn), " not in data directory. Available files: ")
      cat(paste(' -', datalist[-1], collapse="\n" ))
    }
  }
  for (dn in dnames) loadme(dn)
}
freemem <- function (verbose=FALSE) {
  gc()
  if (verbose) {
    message("Garbage collected.")
    print(gc())
    mem <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE))
    mem.gb <- round(mem/1000000, 1)
    message('Free system memory: ', mem.gb, ' GB')
  }
}
makeDataPaths <- setDataPaths #because I always forget which is which
message("   ...Basic functions loaded.")

# ******************************************************************************
#                             PLOTTING FUNCTIONS
# ******************************************************************************

zoomTo <- function(shape, buffer=2) {
  bb <- st_bbox(shape)
  bb[c('xmin', 'ymin')] <- bb[c('xmin', 'ymin')] - buffer
  bb[c('xmax', 'ymax')] <- bb[c('xmax', 'ymax')] + buffer
  coord_sf(xlim=bb[c('xmin', 'xmax')],
           ylim=bb[c('ymin', 'ymax')])
}

# load protected areas maps
# setDataPaths('khaudum_geography', verbose=FALSE)
# load(here::here(procpath, "geographicData.RData"))
# khauPlot <- function(fill="transparent", lwd=1, ...) {
#   ggplot(data=khau) + 
#     geom_sf(fill=fill, linewidth=lwd, ...) + 
#     plot.theme
# }
# plotting hist and params for step length
plotParams <- function(par1, par2, type="step") {
  createWC <- function(x, k, mu, rho) CircStats::dwrpcauchy(x, mu=mu, rho=rho) * k
  createVM <- function(x, k, mu, kap) circular::dvonmises(x, mu=mu, kappa=kap) * k
  createGamma <- function(x, k, mu, sd) {
    gam <- gammaParamsConvert(mean=mu,sd=sd)
    dgamma(x, shape=gam$shape, scale=gam$scale) * k
  }
  if (type == "gamma") {
    paramFun <- createGamma
    data <- pdata$step
  } else if (type == "vonmises") {
    paramFun <- createVM
    data <- pdata$angle
  } else if (type == "wrpcauchy") {
    paramFun <- createWC
    data <- pdata$angle
  } else {
    message("WARNING: 'type' should be one of: gamma, vonmises, wrpcauchy")
    stop()
  }
  data <- data[!is.na(data)]
  h = hist(data, breaks=200)
  k = diff(h$mids[1:2]) * length(data)
  x = seq(min(data), max(data),length=400)
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
           "#0072B2", "#D55E00", "#CC79A7")
  
  for (i in 1:length(par1) ) {
    lines(x, paramFun(x, k, par1[i], par2[i]), col=pal[i], lwd=2)
  }
  legend(max(data)*0.5, max(h$counts)*0.8, 
         legend=c("resting", "foraging", "correlated walk"),
         col=pal, lty=1, cex=0.8, bty="n")
}

plotBase <- function(path, activity=NULL, seed=1001, b=5000) {
  if (is.null(activity)) {
    base <- ggplot()
    title <- paste0("Simulated ele path, seed = ", seed)
  } else {
    color.ref <- list(forage="Greens", move="Reds", shelter="Purples")
    if (activity == "move") {rast = move.rast
    } else if (activity == "shelter") {rast = shelt.rast
    } else if (activity == "forage") {rast = forage.rast
    }
    # subset by extent
    ext = terra::ext(c(min(path$x)-b, max(path$x)+b, 
                       min(path$y)-b, max(path$y)+b))
    rast.sub <- terra::crop(rast, ext)
    
    # turn into data frame and plot
    rast.df = terra::as.data.frame(rast.sub, xy=TRUE)
    names(rast.df) <- c('x', 'y', 'value')
    title = paste0(toupper(activity), ", seed = ", seed)
    base <- ggplot() +
      geom_raster(data = rast.df,
                  mapping=aes(x=x, y=y, fill=value),
                  alpha=0.5) + 
      scale_fill_distiller(palette=color.ref[[activity]], direction=1, name=activity)
  }
  p <- base + 
    theme(axis.title=element_blank(),
          axis.text=element_text(size=18)) +
    ggtitle(title)
  return(p)
}
plotPaths <- function(simRes, barriers, perm=NULL, seed=1001, b=5000, ping="1 hour",
                      activity=NULL, inx=NULL, colorby=NULL) {
  ## debugging
  # seed=1001
  # b=2000
  # activity='move'
  # inx=NULL
  # colorby="inx"
  
  # get simulation locations
  path <- simRes$locations
  if (!is.null(inx)) path <- path[inx,]
  path$INX <- 1:nrow(path)
  xlim = c(min(path$x)-b, max(path$x)+b)
  ylim = c(min(path$y)-b, max(path$y)+b)
  
  # subset path to match GPS data collection methods
  # working under the assumption that i == minute
  if (ping == "1 hour") pings <- seq(1, nrow(path), by=60)
  else if (ping == "5 hours") pings <- seq(1, nrow(path), by=300)
  else if (is.null(ping)) pings = 1:nrow(path)
  path <- path[pings,]
  
  # get home range data
  inputs <- simRes$inputs
  # homerange <- data.frame(x=inputs$in_home_x,
  #                         y=inputs$in_home_y)
  # home_xy = st_as_sf(homerange, coords=c('x', 'y'), crs=32735)
  # home_circ <- st_buffer(home_xy, dist=inputs$in_home_r)
  
  
  # Set up barrier data, adding permeability if applicable
  b.df <- bind_rows(barriers) %>% as.data.frame()
  if (!is.null(perm)) {
    for (i in 1:length(barriers)) barriers[[i]]$perm <- perm[i]
    b.df <- b.df %>% mutate(perm=factor((1-perm)*100))
  }
  b.df <- st_as_sf(b.df)
  
  ###### PLOTTING! ######
  # STEP 1: plot barriers, homerange, and starting point on raster
  pb <- plotBase(path, activity, seed=seed, b=b)
    # geom_point(data=path[1,], mapping=aes(x=x, y=y), color="red", alpha=0.3, size=10) +
    # geom_sf(data=home_circ, alpha=0.2, color='gray', linewidth=2)
  if (!is.null(perm)) {
    p1 <- pb +
      geom_sf(data=b.df, linewidth=2, mapping=aes(alpha=perm)) + 
      guides(alpha=guide_legend("% movement \nblocked"))
  } else {p1 <- pb + geom_sf(data=b.df, linewidth=2)}
  p1 <- p1 + coord_sf(xlim=xlim, ylim=ylim)
  
  # STEP 2: plot paths based on attributes (index, activity, none)
  if (is.null(colorby)) {
    p2 <- p1 + geom_path(data=path, mapping=aes(x=x, y=y), color='black', linewidth=1)
  } else if (colorby == "inx") {
    p2 <- p1 + 
      geom_path(data=path, mapping=aes(x=x, y=y, color=INX), linewidth=1, alpha=0.75) +
      scale_color_distiller(palette='Spectral', direction=1)
  } else if (colorby == "activity") {
    path$behave <- factor(path$behave)
    p2 <- p1 + 
      geom_path(data=path, mapping=aes(x=x, y=y, color=behave), linewidth=1, alpha=0.75) +
      scale_color_brewer(palette='Dark2', direction=1)
  }
  
  # STEP 3: add on starting point, foraging destinations, and shelter locations
  dests <- path %>%
    mutate(grp = paste(destination_x, destination_y)) %>%
    group_by(grp, destination_x, destination_y) %>%
    summarize(.groups="drop") %>% dplyr::select(-grp) %>%
    rename(x=destination_x, y=destination_y)
  p3 <- p2 + 
    geom_point(data=path[1,], aes(x=x, y=y), color='black', size=5) +
    geom_point(data=path[1,], aes(x=x, y=y), color='white', size=3)
    # geom_point(data=dests, aes(x=x, y=y), color="black", size=0.5, shape=15) +
    # geom_point(data=ELE_shelterLocs, aes(x=x, y=y), color="purple", size=4)
  p3
}

message("   ...Plotting functions loaded.")


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                            FUNCTIONS FOR SIMULATIONS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

runSim <- function( barriers, perm, seed=1001, timesteps=5000,
                    activities=NULL, colorby=NULL, thin=0.5 ) {
  set.seed(seed)
  message('Starting simulation.\nseed ', seed)
  simRes <- abmFences::abm_simulate(
    start = start, 
    timesteps = timesteps, 
    des_options = des_options,
    options = options,
    shelterLocations = ELE_shelterLocs,
    shelterSize = ELE_shelterSize,
    avoidPoints = ELE_avoidLocs,
    destinationRange = ELE_destinationRange,
    destinationDirection = ELE_destinationDirection,
    destinationTransformation = ELE_destinationTransformation,
    destinationModifier = ELE_destinationModifier, 
    avoidTransformation = ELE_avoidTransformation,
    avoidModifier = ELE_avoidModifier,
    k_step = ELE_k_step,
    s_step = ELE_s_step,
    mu_angle = ELE_mu_angle,
    k_angle = ELE_k_angle, 
    rescale_step2cell = ELE_rescale,
    behave_Tmat = ELE_behaveMatrix,
    rest_Cycle = ELE_rest_Cycle,
    additional_Cycles = ELE_additional_Cycles,
    shelteringMatrix = ELE_shelter,
    foragingMatrix = ELE_forage,
    movementMatrix = ELE_move,
    barrier=ELE_barriers
  )
  message('simulation complete.')
  
  message('thinning samples...')
  NR = nrow(simRes$locations)
  inx = sample(1:NR, size=thin*NR)
  
  if (activities == "all") activities <- c('move', 'forage', 'shelter')
  if (is.null(activities)) { 
    p <- plotPaths(simRes, barriers, perm, seed=seed, colorby=colorby, inx=inx)
  } else {
    plot_list <- list(
      move   = plotPaths(simRes, barriers, perm, seed, activity = 'move',   colorby=colorby, inx=inx),
      forage = plotPaths(simRes, barriers, perm, seed, activity = 'forage', colorby=colorby, inx=inx),
      shelter= plotPaths(simRes, barriers, perm, seed, activity = 'shelter',colorby=colorby, inx=inx) 
    )
    p <- ggpubr::ggarrange(plotlist = plot_list[activities], nrow=1, common.legend=TRUE)
  }
  print(p)
}

message("   ...Simulation functions loaded.")


# ******************************************************************************
#                          OTHER CUSTOM FUNCTIONS
# ******************************************************************************

projectMe <- function(obj, crs) {
  mycrs = terra::crs(obj, describe=T)$code
  if (mycrs != crs) {
    if (class(obj)[1] == "SpatRaster") obj = terra::project(obj, paste0("EPSG:", crs)) 
    if (class(obj)[1] == "sf") obj = sf::st_transform(obj, crs) 
  }
  return(obj)
}

cleanUpModel <- function(fit) {
  # a function to clean up all the bloat in model returns.
  message("Cleaning up your model! Starting size is: ",
          format(object.size(fit), unit="Mb"))
  message(' removing data...')
  fit$data <- NULL
  fit$y <- NULL
  fit$linear.predictors <- NULL
  fit$weights <- NULL
  fit$fitted.values <- NULL
  fit$model <- NULL
  fit$prior.weights <- NULL
  fit$residuals <- NULL
  fit$effects <- NULL
  fit$qr$qr <- NULL
  message('Ending size is: ', format(object.size(fit), unit="Mb"))
  gc()
  return(fit)
}
message("   ...Custom functions loaded.")


# ******************************************************************************
message("All utility functions and parameters loaded.")
message("Use setDataPaths() to set paths to a specific data directory")
# ******************************************************************************
