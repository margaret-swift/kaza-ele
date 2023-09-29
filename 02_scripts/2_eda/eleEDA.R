# eleEDA.R
# Created 26 Sept 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/eleEDA.R')
source('02_scripts/utilities.R')
setDataPaths('geographic')
load(here(procpath, 'geographic.rdata'))
setDataPaths('elephant')
load(here(procpath, 'ele.rdata'))

# ******************************************************************************
#                                     STATS
# ******************************************************************************

ele.df %>% 
  group_by(ID, SEX) %>% 
  mutate(MYSPEED = abs(DIST) / abs(as.numeric(DTM))) %>% 
  dplyr::select(SPEED)

# ******************************************************************************
#                                    PLOTTING
# ******************************************************************************

plotPath <- function(i) {
  data <- ele.lines[i,]
  title = paste0('Elephant ', data$ID, ' (', data$SEX, ')')
  bb <- st_bbox(data)
  cols <- ifelse(data$SEX == "F", 'black', 'orange')
  ggplot() +
    geom_sf(data=waters_nat, color='lightblue') + 
    geom_sf(data=waters_art, color='darkblue', size=3) + 
    geom_sf(data=data, mapping=aes(color=SEX)) +
    geom_sf(data=fences, fill=NA, color='red', linewidth=2) +
    ggtitle(title) + plot.theme +
    scale_color_manual(values=cols) +
    guides(color="none") +
    coord_sf(xlim=c(bb['xmin'], bb['xmax']), 
             ylim=c(bb['ymin'], bb['ymax']))
}
makeHist <- function(i) {
  data <- ele.df %>% 
    filter(ID == i) %>% 
    mutate(DATE = date(DATE.TIME))
  hist.data <- data %>% 
    group_by(DATE) %>% 
    summarize(n=n()) %>%
    mutate(FLAG = ifelse(n > 27, "HIGH", ifelse(n<20, "LOW", "AVG")))
  title = paste0('Elephant ', data$ID[1], ' (', data$SEX[1], ')')
  cols <- list(AVG='darkgray', HIGH='#08c952', LOW='#f2055c')[unique(hist.data$FLAG)]
  ggplot() + 
    geom_bar(data=hist.data, 
             mapping=aes(x=DATE, y=n, fill=FLAG), 
             stat="identity") + 
    ggtitle(title) + 
    scale_fill_manual(values=cols) + 
    plot.theme + 
    ylab('# fixes per day')
}

# Plots for all IDs
ids <- unique(ele.df$ID)
for (i in ids) {
  ggsave(filename = here(outdir, "eda", "ele_freq",
                         paste0("ele_", i, ".png")),
         plot=makeHist(i),
         width=10, height=5)
  ggsave(filename = here(outdir, "eda", "ele_paths", 
                         paste0("ele_", i, ".png")),
         plot=plotPath(i))
}


# plot all histograms together that have too many points
data <- ele.df %>% 
  mutate(DATE = date(DATE.TIME)) %>% 
  group_by(ANIMAL_ID, DATE) %>% 
  summarize(n=n()) %>%
  mutate(FLAG = ifelse(n > 27, "HIGH", ifelse(n<20, "LOW", "AVG")))
cols <- list(AVG='darkgray', HIGH='#08c952', LOW='#f2055c')

my_ids <- data %>% 
  ungroup() %>% 
  group_by(ANIMAL_ID, FLAG) %>% 
  summarize(n=n()) %>% 
  filter(FLAG == "HIGH") %>% 
  pull(ANIMAL_ID)
  
ggplot() + 
  geom_bar(data=data %>% filter(ANIMAL_ID %in% my_ids), 
           mapping=aes(x=DATE, y=n, fill=FLAG), 
           stat="identity") + 
  scale_fill_manual(values=cols) + 
  facet_wrap(~ANIMAL_ID, scales="free_y", ncol=1) +
  plot.theme + 
  ylab('# fixes per day')


makeHistAll()

