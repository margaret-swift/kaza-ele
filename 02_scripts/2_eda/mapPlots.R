# mapPlots.R
# Created 10 January 2024
# Margaret Swift <mes114@duke.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

# tell R where you are & load utilities
here::i_am('02_scripts/2_eda/mapPlots.R')
source(here::here('02_scripts','utilities.R')) #this is for other custom funcs

# load the "world" dataset from spData
library(spData)

# pull geographic data (KAZA boundary, fences, roads, etc.)
quickload()
aoi_utm <- st_transform(kaza_aoi, crs="EPSG:32734")
aoi_alb <- st_transform(kaza_aoi, crs="ESRI:102022")


# where to save
outfile = "./03_output/maps"


# ******************************************************************************
#                                 FUNCTIONS & PARAMETERS
# ******************************************************************************

# function to zoom to a shape with a user-provided buffer
zoomTo <- function(shape, buffer=NULL, buffer.x=NULL, buffer.y=NULL) {
  if (is.null(buffer.x)) buffer.x <- buffer
  if (is.null(buffer.y)) buffer.y <- buffer
  bb <- st_bbox(shape)
  bb[c('xmin', 'ymin')] <- bb[c('xmin', 'ymin')] - buffer.x
  bb[c('xmax', 'ymax')] <- bb[c('xmax', 'ymax')] + buffer.y
  coord_sf(xlim=bb[c('xmin', 'xmax')],
           ylim=bb[c('ymin', 'ymax')])
}

# function to add transparent background to plot
#   NOTE: figure must be saved as PNG to have transparent BG!
transparentBg <- function() {
  theme(
    axis.title.x=element_blank(),
    axis.text.x =element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y =element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.background= element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )
}

## GEOGRAPHIC OBJECTS
sadc_list <- c("south africa", "lesotho", "eswatini", "angola",
               'mozambique', 'malawi', 'zambia', 'zimbabwe', 
               'botswana', 'namibia', 'tanzania', 'mauritius', 
               'democratic republic of the congo' , 'madagascar')
kaza_members <- c('angola', 'zambia', 'zimbabwe', 'botswana', 'namibia')
africa <- world %>% 
  filter(continent == "Africa", !is.na(iso_a2)) %>% 
  st_transform(crs=st_crs(fences))
sadc <- africa %>% filter(tolower(name_long) %in% sadc_list)
sadc_kaza <- africa %>% filter(tolower(name_long) %in% kaza_members)
bots <- africa %>% filter(name_long=="Botswana")

# ******************************************************************************
#                                    PLOTTING
# ******************************************************************************

## ROADS
ggplot() +
  geom_sf(data=roads, color='black', linewidth=2) +
  geom_sf(data=roads, color='white', linewidth=1) +
  transparentBg() +
  zoomTo(kaza, buffer=1)
ggsave(filename=file.path(outfile, 'roads.png'),
       width=10, height=10, bg='transparent')

## RIVERS
ggplot() +
  geom_sf(data=rivers,color='blue', linewidth=1) +
  transparentBg() +
  zoomTo(kaza, buffer=1)
ggsave(filename=file.path(outfile, 'rivers.png'),
       width=10, height=10, bg='transparent')

## FENCES
ggplot() +
  geom_sf(data=fences,color='black', linewidth=1) +
  transparentBg() +
  zoomTo(kaza, buffer=1)
ggsave(filename=file.path(outfile, 'fences.png'),
       width=10, height=10, bg='transparent')

## ROADS AND FENCES ON TOP OF KAZA MAP
ggplot() +
  geom_sf(data=kaza, fill='green', alpha=0.5) +
  geom_sf(data=fences, color='black', linewidth=1) +
  geom_sf(data=roads, color='black', linewidth=2) +
  geom_sf(data=roads, color='white', linewidth=1) +
  transparentBg() +
  zoomTo(kaza, buffer=1)
ggsave(filename=file.path(outfile, 'linear_features_kaza.png'),
       width=10, height=10, bg='transparent')

## SADC REGION WITH KAZA
ggplot() + 
  geom_sf(data=sadc) + 
  geom_sf(data=kaza, fill='green', alpha=0.5, color=NA) +
  transparentBg()
ggsave(filename=file.path(outfile, 'sadc_kaza.png'),
       width=10, height=10, bg='transparent')

## KAZA WITH MEMBER COUNTRIES
ggplot() + 
  geom_sf(data=africa, fill="darkgray", color="#868686") + 
  geom_sf(data=sadc_kaza, fill='#e6e6e6') + 
  geom_sf(data=kaza, fill='green', alpha=0.5, color=NA) +
  transparentBg() + 
  zoomTo(kaza, buffer=7.5)
ggsave(filename=file.path(outfile, 'kaza_members.png'),
       width=10, height=10, bg='transparent')

## AFRICA PLAIN (LIGHT)
ggplot() + 
  geom_sf(data=africa, fill='#e6e6e6', color=NA) + 
  transparentBg()
ggsave(filename=file.path(outfile, 'africa_light.png'),
       width=10, height=10, bg='transparent')

## AFRICA PLAIN (DARK)
ggplot() + 
  geom_sf(data=africa, fill='darkgray', color=NA) + 
  transparentBg()
ggsave(filename=file.path(outfile, 'africa_dark.png'),
       width=10, height=10, bg='transparent')

## KHAUDUM WITH RIVERS
ggplot() + 
  geom_sf(data=khau) + 
  geom_sf(data=rivers %>% filter(BB_DIS_ORD<=6), color='blue', linewidth=2) +
  transparentBg() + 
  zoomTo(khau, buffer=0)
ggsave(filename=file.path(outfile, 'khau.png'),
       width=10, height=10, bg='transparent')


## AREA OF INTEREST
library(tmaptools)
# load(file=here(outdir, "hmm", "hmm_ssf.rdata"))
# bb <- st_bbox(data) %>% bb_poly()
# st_write(bb, dsn=rawpath('namib_ele_aoi', "namib_ele_aoi.shp"))

# data = hmm.ssf.dat %>% st_as_sf(coords=c('x', 'y'), crs=32734)
# ggplot() + 
  # geom_sf(data=aoi, fill=NA, color='black', linewidth=1) + 
  # geom_sf(data=bb, fill=NA, color='green')
  # geom_sf(data=data %>% sample_n(50000), size=0.2, alpha=0.2)



# ******************************************************************************
#                               PLOTTING RASTERS
# ******************************************************************************


## METADATA
setDataPaths('landcover')
lm05 <- read.csv(metapath('landcover_meta_2005.csv'), )
struc.col <- data.frame(
  value = 1:8,
  structure = unique(lm05$structure),
  colors=c('lightgray',
           '#17b39c',
           '#6de3d2',
           'darkgreen', 
           "purple",
           "#B5DC83",
           "#66A61E",
           "#E6AB02")
  )


## WWF LANDCOVER TYPES -- 2021 classification
##  https://storymaps.arcgis.com/stories/db69e30bc8434330ab9445a200b19709
# lm21 <- read.csv(metapath('landcover_meta_2021.csv'), )
# lr21 <- terra::rast(rawpath('eoss4wwf_kaza_tfca_landcover_2020.tif'))
# names(lr21) <- 'category'
# lr21.crop <- terra::crop(lr21, aoi_utm)
# rc21 <- lm21 %>% 
#   left_join(struc.col, by="structure") %>% 
#   dplyr::select(category, value)
# lr21.re <- terra::classify(lr21.crop, rc21)
# terra::plot(lr21.re,
#             col=struc.col$colors)


## WWF LANDCOVER TYPES -- 2005 classification
##  https://www.arcgis.com/home/item.html?id=b82b036653df4ef8abf184f1c4c656ec

lr05 <- terra::rast(rawpath('landcover_2005_fixed.tif'))
lr05.crop <- terra::crop(lr05, aoi_alb)
rc05 <- lm05 %>% 
  left_join(struc.col, by="structure") %>% 
  dplyr::select(category, value)
lr05.re <- terra::classify(lr05.crop, rc05)
cols <- struc.col[struc.col$value %in% unique(lr05.re)$agg_class,c(1,3)]
terra::plot(lr05.re, col=cols)


## EVI
setDataPaths('evi')
library(terra)

erasts <- list.files(rawpath('evi_kaza_aoi'), pattern="2019", full.names=TRUE)
er.dry.1 <- rast(erasts[1])
er.dry.2 <- rast(erasts[2])
er.wet.1 <- rast(erasts[3])
er.wet.2 <- rast(erasts[4])

min = -0.1; max = 0.30; ncol=8
vals <- seq(min, max, length.out=ncol+1)
color.ramp <- data.frame(
  from = vals[1:(ncol)],
  to = vals[2:(ncol+1)],
  color = RColorBrewer::brewer.pal(ncol, 'YlGn')
)

terra::plot(er.wet.1, type="continuous", col=color.ramp, 
            main="wet season 1, 2016 (Oct-Jan 2015)")
terra::plot(er.wet.2, type="continuous", col=color.ramp, 
            main="wet season 2, 2016 (Jan-May 2016)")
terra::plot(er.dry.1, type="continuous", col=color.ramp, 
            main="dry season 1, 2016 (May-July 2016)")
terra::plot(er.dry.2, type="continuous", col=color.ramp, 
            main="dry season 2, 2016 (July-Oct 2016)")
