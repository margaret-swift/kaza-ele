# utilities.R
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                    LIBRARY LOADING & DEFINE GLOBAL OBJECTS
# ******************************************************************************
message("\nLoading all utility functions and parameters from utilities.R\n")

message('Loading base packages for all projects...')
pacman::p_load(tidyverse, patchwork, 
               sf, rgdal, 
               lubridate, here)
message('   ...Packages loaded.')

message("Setting base objects...")
datadir <- here::here("01_data")
scriptdir <- here::here("02_scripts")
outdir <- here::here("03_output")
crs <- CRS("+proj=longlat +zone=34S +datum=WGS84")
message("   ...Base objects loaded.")

# Plot Params
message("Setting plot params...")
plot.theme <- theme( text = element_text(size=18), 
                     title = element_text(size=20),
                     axis.text.x = element_text(angle=45, vjust=0.5))
blank.theme <- theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
scalebar <- ggspatial::annotation_scale( pad_x = unit(0.05, "in"), 
                                         pad_y = unit(0.05, "in"),
                                         text_face="bold") 
message("   ...Plot params loaded.")

# ******************************************************************************
#                             BASE UTILITY FUNCTIONS
# ******************************************************************************
message("Loading functions...")

as.c <- as.character
as.n <- as.numeric
as.v <- as.vector
nog <- st_drop_geometry
len <- length

removePunct <- function(v) gsub('\\)|\\(|\\/|\\.|\\?', '', v)
replaceSpace <- function(v, r="_") gsub(' ', r, v)
fixNames <- function(d) {
  names(d) <- replaceSpace(removePunct(toupper(names(d))))
  d
}

makeDataPaths <- function(dataname, verbose=TRUE, print.length=5) {
  datapath <<- file.path(datadir, dataname)
  files <- list.files(datapath)
  if (!length(files)) {
    warning('No files found in data directory "', dataname, 
            '".\nPossible data file names: ',
            paste("\n     -", list.files(datadir)))
  } else {
    metapath <<- file.path(datapath, 'meta')
    rawpath  <<- file.path(datapath, 'raw')
    procpath <<- file.path(datapath, 'processed')
    
    if (verbose) {
      printPathFiles <- function(path) {
        cat(gsub(".*01_data", "", path))
        files <- list.files(path)
        L=min(length(files), print.length)
        message(paste("\n     -", files[1:L]))
        if (length(files)>print.length) message('     ...')
      }
      cat('resetting data paths to...\n')
      printPathFiles(metapath)
      printPathFiles(rawpath)
      printPathFiles(procpath)
    }
  }
}
setDataPaths <- makeDataPaths
message("   ...Basic functions loaded.")

# ******************************************************************************
#                             PLOTTING FUNCTIONS
# ******************************************************************************

# load protected areas maps
# setDataPaths('khaudum_geography', verbose=FALSE)
# load(here::here(procpath, "geographicData.RData"))
# khauPlot <- function(fill="transparent", lwd=1, ...) {
#   ggplot(data=khau) + 
#     geom_sf(fill=fill, linewidth=lwd, ...) + 
#     plot.theme
# }

message("   ...Plotting functions loaded.")


# ******************************************************************************
#                          OTHER CUSTOM FUNCTIONS
# ******************************************************************************

message("   ...Custom functions loaded.")


# ******************************************************************************
#                          ADVANCED HOUSEKEEPING
# ******************************************************************************

# Set of functions to list all custom functions and variables
listFuns <- function() {
  df <- my.funs
  types <- levels(df$ftype)
  for (type in types) {
    cat(type, "\n", df$fname[df$ftype == type], "\n\n")
  }
}
listVars <- function() {
  df <- my.vars
  types <- levels(df$vtype)
  for (type in types) {
    cat(type, "\n", df$vname[df$vtype == type], "\n\n")
  }
}
.listVars <- function() {
  lst <- ls(envir = .GlobalEnv)
  inx <- sapply(lst,function(var) any(class(get(var))!='function'))
  vs <- lst[inx]
  
  getType <- function(val) {
    val <- tolower(val)
    type <- "OTHER"
    if (grepl("path", val)) { type = "PATH"
    } else if (grepl("theme|scale|crs", val)) { type = "THEME"
    } else if (grepl("khau|kaza|waters", val)) { type = "MAPDATA"
    } 
    type = factor(type, levels=c("PATH", "THEME", "MAPDATA", "OTHER"))
  }
  vtypes <- unlist(lapply(vs, getType))
  df <- data.frame(vtype=vtypes, vname=vs) %>% arrange(vtype)
  df
}
.listFuns <- function() {
  lst <- ls(envir = .GlobalEnv)
  inx <- sapply(lst,function(var) any(class(get(var))=='function'))
  fs <- lst[inx]
  
  getType <- function(val) {
    val <- tolower(val)
    type <- "OTHER"
    if (grepl("as.", val)) { type = "LAZY"
    } else if (grepl("plot", val)) { type = "PLOTTING"
    } else if (grepl("list", val)) { type = "LISTING"
    } else if (grepl("fix|^re", val)) { type = "FIXING"
    } else if (grepl("set", val)) { type = "HOUSEKEEPING"
    } 
    type = factor(type, levels=c("HOUSEKEEPING", "FIXING", "LISTING", 
                                 "PLOTTING", "LAZY", "OTHER"))
  }
  ftypes <- unlist(lapply(fs, getType))
  df <- data.frame(ftype=ftypes, fname=fs) %>% arrange(ftype)
  df
}

# only create list of custom vars and funs once, via housekeeping switch
if ( !exists('housekeeping_switch') ) {
  my.vars <- .listVars()
  my.funs <- .listFuns()
  housekeeping_switch = TRUE
} 

message("   ...Advanced housekeeping functions loaded.\n")

# ******************************************************************************
message("All utility functions and parameters loaded.")
message("Use listFuns() to list all available utility functions, 
        and listVars() to list all base parameters.")
# ******************************************************************************
