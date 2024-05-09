# elePlots.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/2_eda/elePlots.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(sp, adehabitatHR, reshape2, 
               spData, ggspatial, move, terra)
library(moveVis)
setDataPaths('geographic')
load(procpath('geographic.rdata'))
er <- rast(rawpath('khaudum_landsat_rasters', 
                   'waterdist_evi_raster_dry_2017.tif'))
setDataPaths('elephant')
load(procpath('ele.rdata'))


# ******************************************************************************
#                                    PLOTTING
# ******************************************************************************
makeHist <- function(i) {
  data <- ele.df %>% nog() %>% filter(ID == i)
  hist.data <- data %>%
    group_by(DATE, .drop=FALSE) %>%
    summarize(n=n()) %>%
    mutate(FLAG = ifelse(n > 27, "HIGH", ifelse(n<20, "LOW", "AVG")))
  title = paste0('Elephant ', data$ID[1], ' (', data$SEX[1], ')')
  cols <- list(AVG='darkgray', HIGH='#08c952', LOW='#f2055c')[unique(hist.data$FLAG)]
  ggplot() +
    geom_bar(data=hist.data,
             mapping=aes(x=DATE, y=n, fill=FLAG),
             stat="identity") +
    ggtitle(title) + theme(axis.title.x=element_blank()) +
    scale_fill_manual(values=cols) + ylab('number of fixes') +
    scale_x_date(breaks='1 year', date_labels = "%Y") + theme(text=element_text(size=20))
}
plotXY <- function(i) {
  data <- ele.df %>% nog() %>% 
    filter(ID == i) %>% 
    mutate(
      BURST = factor(BURST),
      MONTH = month(DATE.TIME),
      YEAR = year(DATE.TIME),
      MONTHYEAR = paste(MONTH, YEAR),
      SEASON = ifelse(MONTH %in% 4:10, "DRY", "WET"))
  pbase <- ggplot(data=data, 
                  mapping=aes(x=as.Date(DATE.TIME), 
                              group=MONTHYEAR, color=SEASON)) +
    plot.theme + 
    xlab('') + 
    scale_color_brewer(palette="Dark2", direction=-1)+
    scale_x_date(date_labels = "%m-%Y", breaks='3 months')
  
  px <- pbase +
    geom_line(mapping=aes(y=X)) #+
  geom_point(data=data %>% filter(is.na(FIXRATE)),
             aes(X), color='black', size=2) +
    ggtitle(paste(i, data$SEX[1]))
  py <- pbase + geom_line(mapping=aes(y=Y)) +
    geom_point(data=data %>% filter(is.na(FIXRATE)),
               aes(y=Y), color='black', size=2)
  px / py
}
plotRegion <- function(i, t=z) {
  bb <- ele.df %>% filter(ID==i) %>% st_bbox()
  tt <- 0.1
  ggplot() + 
    geom_sf(data=kaza, linewidth=1, fill='#94EF81') +
    geom_sf(data=fences, color='black', linewidth=1) +
    geom_sf(data=rivers, color='blue', linewidth=1) +
    geom_sf(data=roads, color='black', linewidth=2) + 
    geom_sf(data=roads, color='white', linewidth=1) + 
    # annotate(geom = "rect", 
    #          xmin=bb['xmin']-tt, xmax=bb['xmax']+tt,
    #          ymin=bb['ymin']-tt, ymax=bb['ymax']+tt,
    #          colour = "red", fill = NA, linewidth=1) +
    coord_sf(xlim=c(bb['xmin']-t, bb['xmax']+t),
             ylim=c(bb['ymin']-t, bb['ymax']+t)) + 
    theme(text=element_text(size=20)) + 
    annotation_scale(text_cex=1)
}
makeY <- function(i) {
  data <- ele.df %>%
    group_by(ID, SEX) %>%
    filter(ID == i, as.n(DTM) > 45) %>%
    mutate(
      BURST = factor(BURST),
      MONTH = month(DATE.TIME),
      YEAR = year(DATE.TIME),
      MONTHYEAR = paste(MONTH, YEAR),
      MPH = DIST / (as.n(DTM)/60))
  base = ggplot( data=data,
          mapping=aes(x=as.Date(DATE.TIME), group=MONTHYEAR, color=SEASON)) +
    theme(axis.title.x=element_blank()) +
    scale_color_brewer(palette="Dark2", direction=-1) +
    scale_x_date(breaks='1 year', date_labels = "%Y") + theme(text=element_text(size=20))
  base
}
plotPoints <- function(i) {
  i=20
  data <- ele.df %>% filter(ID==i) %>% sample_frac(0.5)
  # data <- data %>% filter(DATE > as.Date("2021-01-01"), 
                          # DATE < as.Date("2021-02-28"))
  bb <- st_bbox(data)
  width=abs(bb[['xmin']] - bb[['xmax']])
  height=abs(bb[['ymin']] - bb[['ymax']])
  nrow = ifelse(width > height, 2, 1)
  
  ext <- bb %>% 
    st_as_sfc() %>% 
    st_transform("EPSG:32734") 
  ggplot() + 
    geom_sf(data=data, alpha=0.1, mapping=aes(color=DATE)) +
    geom_sf(data=fences, color='black') + 
    geom_sf(data=roads, color='black', linewidth=2) + 
    geom_sf(data=roads, color='white', linewidth=1) + 
    geom_sf(data=rivers, color='blue', linewidth=2) + 
    geom_sf(data=waters_art, color='darkblue', size=3) +
    facet_wrap(~SEASON, nrow=nrow) +
    coord_sf(xlim=c(bb['xmin'], bb['xmax']),
             ylim=c(bb['ymin'], bb['ymax'])) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) + 
    annotation_scale(text_cex=1) + 
    plot.theme + blank.theme
}
plotAll <- function(i) {
  p1 <- makeHist(i) 
  p2 <- makeY(i) + geom_line(mapping=aes(y=X))
  p3 <- makeY(i) + geom_line(mapping=aes(y=Y))
  p4 <- plotRegion(i)
  p5 <- plotPoints(i) 
  
  ( ( p1 + p2 + p3 + plot_layout(ncol=1) ) | 
      ( p4 + p5 + plot_layout(ncol=1, heights=c(1,2)) ) ) + 
    plot_layout( widths=c(3,1))
}

plotPoints(5)

# ******************************************************************************

# Plot all and save
ids <- unique(ele.df$ID)
# tic()
# for (i in ids) {
#   message(i)
#   p=plotAll(i)
#   ggsave(filename = here(outdir, "eda", "ele_all",
#                          paste0("ele_", i, ".png")),
#          plot=p, width=21, height=13)
# }
# toc()


# Plot tracks
library(amt)
i=20
plotPoints(i)
plotRegion(i, t=0.5)

# video
move_data <- ele.df %>% 
  mutate(time = as.POSIXct(DATE.TIME),
         YEAR = year(DATE.TIME)) %>% 
  filter(ID==i, 
         DATE > as.Date("2021-02-01"), 
         DATE < as.Date("2021-02-28")) %>% 
  df2move( proj="EPSG:4326",
           x="LON", y="LAT", 
           time="time", 
           track_id="ID"
          ) %>% 
  align_move(res=1, unit="hours")

# create spatial frames with a OpenStreetMap watercolour map
frames <- frames_spatial(move_data, 
                         map_service = "carto", 
                         map_type="light",
                         alpha = 0.5) %>% 
  add_labels(x = "Longitude", y = "Latitude") %>% # add some customizations, such as axis labels
  add_northarrow() %>% 
  add_scalebar() %>% 
  add_timestamps(type = "label") %>% 
  add_progress()

# animate frames
fname = here("03_output", "movies", "ele20_FEB2020.gif")
animate_frames(frames, out_file=fname)


# plot fence crossing
df <- data.frame(sex=c('female', 'male'),
                 river = c(10.1, 14.5),
                 road = c(15.3, 25.8),
                 fence = c(0, 3.5)
                 ) %>% 
  reshape2::melt(id.var='sex')

ggplot(data=df, aes(x=variable, y=value*0.01, fill=sex)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c('black', 'darkgray')) +
  theme(text=element_text(size=24)) + 
  xlab('linear feature') + ylab('chance of crossing')






# ******************************************************************************
#                               PLOTTING HMM OUTPUTS
# ******************************************************************************

ids <- unique(ele.khau$ID)
data <- ele.khau %>% 
  mutate(YEAR = year(DATE),
         is_new = STATE != lag(STATE)) %>% 
  filter(ID == 2, 
         DATE > as.Date("2015-01-01"), 
         DATE < as.Date("2015-12-31")) %>% 
  st_as_sf(coords=c("LON", "LAT"), crs=4326) %>% 
  ungroup()
data$NEWBURST <- cumsum( data$is_new )
inx.keep <- which(rle(data$NEWBURST)$lengths > 1)
lines <- data %>% 
  group_by(BURST) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING")
lines_states <- data %>% 
  filter(NEWBURST %in% inx.keep) %>% 
  group_by(NEWBURST, STATE) %>% 
  dplyr::summarize(do_union=FALSE) %>% 
  st_cast("LINESTRING")

ggplot() + 
  geom_sf(data=lines) +
  # geom_sf(data=lines_states, mapping=aes(color=STATE)) +
  # geom_sf(data=data, size=2, alpha=0.75, pch=1) +
  geom_sf(data=data, size=3, mapping=aes(color=STATE)) +
  # geom_sf(data=waters_art) +
  annotation_scale(text_cex=1) +
  theme(text=element_text(size=18)) + 
  guides(color="none") + 
  scale_color_manual(values=c('orange', 'purple', 'green'))

ggplot(data) + 
  geom_histogram(aes(fill=STATE, x=DIST)) + 
  facet_wrap(~STATE, ncol=1, scales="free_y") +
  scale_fill_manual(values=c('orange', 'purple', 'green')) + 
  theme(text=element_text(size=18)) + 
  ylab('frequency') + xlab('step length (m)')

ggplot(data) + 
  geom_histogram(aes(fill=STATE, x=REL.ANGLE)) + 
  facet_wrap(~STATE, ncol=1, scales="free_y") +
  scale_fill_manual(values=c('orange', 'purple', 'green')) + 
  theme(text=element_text(size=18)) + 
  ylab('frequency') + xlab('turning angle')




# ******************************************************************************

# measuring displacement


# moving window

i=1; j=25;
data <- ele.df %>% filter(ID == i)
data.sp <- data %>%
  nog() %>% 
  dplyr::select(X, Y) %>% 
  SpatialPoints()
kernel.ref <- kernelUD(data.sp, h='href')
max = findmax(kernel.ref)
max.pts <- max %>% 
  st_as_sf(crs=32734) %>% 
  st_set_crs(32734) %>% 
  st_transform(crs=st_crs(data)) %>% 
  mutate(SINK = factor(1:length(max)))

whichminList <- function(df) {
  apply(df, 1, which.min)
}
maxdists <- spDists(data.sp, max)
data$SINK <- factor( whichminList(maxdists) )
data$SINK.DIST <- unlist(
  lapply(1:nrow(data), function(e) maxdists[e,data$SINK[e]])
)

# plot distance from sinks
ggplot() + 
  geom_bar(data=data, stat='identity',
           mapping=aes(x=DATE, y=SINK.DIST, 
                       color=SINK, fill=SINK)) 
  # scale_x_datetime(date_breaks="6 months")


# plot with sinks
ggplot() + 
  geom_sf(data=data, mapping=aes(color=SINK), alpha=0.2) +
  geom_sf(data=max.pts, color='black', size=3) +
  geom_sf(data=max.pts, mapping=aes(color=SINK), size=1.5)

# plot UD kernel
plot(kernel.ref)
points(max, pch=19, col='black', cex=1)
points(max, pch=19, col='white', cex=0.5)

#Wayne Getz