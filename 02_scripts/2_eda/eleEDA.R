# eleEDA.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/2_eda/eleEDA.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(sp, adehabitatHR, reshape2)
setDataPaths('geographic')
load(here(procpath, 'geographic.rdata'))
setDataPaths('elephant')
load(here(procpath, 'ele.rdata'))


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
plotRegion <- function(i) {
  bb <- ele.df %>% filter(ID==i) %>% st_bbox()
  t=2
  ggplot() + 
    geom_sf(data=kaza, linewidth=1, fill=NA) +
    geom_sf(data=fences) +
    geom_sf(data=rivers, color='purple') +
    geom_sf(data=roads, color='black', linewidth=1.5) + 
    geom_sf(data=roads, color='white', linewidth=0.5) + 
    annotate(geom = "rect", 
             xmin=bb['xmin'], xmax=bb['xmax'],
             ymin=bb['ymin'], ymax=bb['ymax'],
             colour = "red", fill = NA, linewidth=1) +
    coord_sf(xlim=c(bb['xmin']-t, bb['xmax']+t),
             ylim=c(bb['ymin']-t, bb['ymax']+t))
}
plotPoints <- function(i) {
  data <- ele.df %>% filter(ID==i) %>% sample_frac(0.5)
  bb <- st_bbox(data)
  width=abs(bb[['xmin']] - bb[['xmax']])
  height=abs(bb[['ymin']] - bb[['ymax']])
  nrow = ifelse(width > height, 2, 1)
  ggplot() + 
    # geom_sf(data=waters_nat %>% sample_frac(0.5), color='lightblue') +
    geom_sf(data=data,
            mapping=aes(color=SEASON), alpha=0.1) +
    geom_sf(data=fences, color='black') + 
    geom_sf(data=roads, color='black', linewidth=1) + 
    geom_sf(data=roads, color='white', linewidth=0.5) + 
    geom_sf(data=rivers, color='purple') + 
    geom_sf(data=waters_art, color='darkblue', size=3) +
    scale_color_brewer(palette="Dark2", direction=-1) + 
    facet_wrap(~SEASON, nrow=nrow) +
    coord_sf(xlim=c(bb['xmin'], bb['xmax']),
             ylim=c(bb['ymin'], bb['ymax'])) +
    guides(color='none') + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
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

plotPoints(1)

# ******************************************************************************

# Plot all and save
ids <- unique(ele.df$ID)
tic()
for (i in ids) {
  message(i)
  ggsave(filename = here(outdir, "eda", "ele_all",
                         paste0("ele_", i, ".png")),
         plot=plotAll(i), width=21, height=13)
}
toc()


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