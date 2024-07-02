# createMovementMatrices.R
# Created 06 May 2024
# Margaret Swift <margaret.swift@cornell.edu>


# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/3_modeling/habitat_analyses/createMovementMatrices.R')
source(here::here('02_scripts','utilities.R'))
load(here(outdir, 'habitat_selection', 'hsf.rdata'))
pacman::p_load(stars, tmap, tmaptools, sf, terra)


# ******************************************************************************
#                    LOAD, RECLASS AND CROP GEOGRAPHIC DATA
# ******************************************************************************

setDataPaths('geographic')
load(procpath('geographic.rdata'))

# landcover types metadata
lands.meta <- read.csv(rawpath('kaza_landcover', 'landcover_metadata.csv'))
khau.new <- st_transform(khau, "EPSG:32734") %>% as("Spatial")
bb <- st_bbox(khau.new)
bb[1:4] <- c(xmin=48e4, ymin=795e4, xmax=49e4, ymax=796e4) 
bb <- bb %>% st_as_sfc()

# EVI and water dists 
er <- stars::read_stars(rawpath('khaudum_landsat_rasters', 'waterdist_evi_raster_dry_2017.tif'))
er.khau <- crop(er$EVI, khau.new)
er.sub <- crop(er.khau, bb)

# reclassification scheme
dict <- lands.meta[,c('category', 'cover')]
labels <- lands.meta[,c('category', 'cover', 'cover_class.x')] %>% 
  rename(cover_class=cover_class.x) %>% 
  distinct()

# landcover data from WWF

## BOOKMARK : I am REALLY struggling with getting R to reclassify a stars_proxy object; 
##    should I just go back to using terra?
lr <- rast(rawpath('kaza_landcover', 'eoss4wwf_kaza_tfca_landcover_2020.tif'))
names(lr) <- 'category'
lr.khau <- terra::crop(lr, khau.new) %>%  
  terra::classify(rcl = dict, right = NA)
levels(lr.khau) <- labels
lr.sub <- crop(lr.khau, bb)

### BROKEN STARS WAY ###
# lr <- stars::read_stars(rawpath('kaza_landcover', 'eoss4wwf_kaza_tfca_landcover_2020.tif'))
# lr.khau <- st_crop(lr, khau.new) 
# lr.khau <- lr.khau %>% 
#   mutate(x = category,
#          x = case_when (x<100 ~ NA_real_,   x>100 & x<150 ~ 2))


# world settlement footprint
wsf <- rast(rawpath('WSF2019', 'wsf_khau_binary_utm.tif'))
names(wsf) <- 'is_settled'
wsf.khau <- crop(wsf$is_settled, khau.new)
wsf.khau <- ifel(wsf.khau == 0, NA, 1)
wsf.sub <- crop(wsf.khau, bb)

# rivers
rivers <- st_transform(rivers, "EPSG:32734")


# ******************************************************************************
#                               COMPUTE AND PLOT RASTERS
# ******************************************************************************
# https://www.wvview.org/os_sa/15b_Raster_Analysis_terra.html

pal <- rev(hcl.colors(7, "ag_GrnYl"))

p1 <- tm_shape(lr.khau) +
  tm_raster(style= "cat",
            palette=pal,
            title="Land Cover")+
  tm_shape(khau.new) + 
    tm_borders(col="red", lwd=3) + 
  tm_shape(bb) + 
    tm_borders(col='green', lwd=2)

p2 <- tm_shape(lr.sub)+
  tm_raster(style= "cat",
            palette = pal,
            title="Land Cover", 
            drop.levels=FALSE)

p3 <- tm_shape(er.sub) +
  tm_raster(style= "pretty",
            title="Enhanced Vegetation Index", 
            palette="Greens",
            drop.levels=FALSE)
tmap_arrange(p1, p2, p3, ncol=3)

# dist <- terra::distance(wsf.khau, rivers)


# ******************************************************************************
#                               CALCULATE MATRICES
# ******************************************************************************

state <- 'foraging'

biOut <- (lc==42 | lc == 43) & (dem2 > 1000)


cfs <- as.data.frame(t(coefficients(hsf)))




cfs$stateresting
cfs$`stateresting:waterdist`

covers <- t(cfs[,c(grepl('cover_class', names(cfs)))]) %>% 
  as.data.frame() %>% 
  rownames_to_column()
names(covers) <- c('var', 'val')
# covers$trueval <- covers$val + cfs$`(Intercept)`
covers$trueval <- exp(covers$val)
covers$var <- gsub('cover_class', '', covers$var)
ggplot() + 
  geom_bar(data=covers,
           mapping=aes(x=var, y=trueval, fill=var),
           stat='identity') + 
  theme(text=element_text(size=24)) + 
  xlab('landcover type') + 
  ylab('increase in land use over bare ground') + 
  scale_fill_manual(values=pal[-1])



states <- t(cfs[,c(grepl('state', names(cfs)))]) %>% 
  as.data.frame() %>% 
  rownames_to_column()
names(states) <- c('var', 'val')
states$trueval <- exp(states$val)
states$var <- gsub('state', '', states$var)
ggplot() + 
  geom_bar(data=states %>% filter(!grepl(':', var)),
           mapping=aes(x=var, y=trueval),
           stat='identity') + 
  theme(text=element_text(size=24)) + 
  xlab('landcover type') + 
  ylab('increase in land use exploring')






















