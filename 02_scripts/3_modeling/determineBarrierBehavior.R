# determineBarrierBehavior.R
# Created 09 November 2023
# Margaret Swift <mes114@duke.edu>

# Xu, W., Dejid, N., Herrmann, V., Sawyer, H. & Middleton, A. D. 
#   Barrier Behaviour Analysis (BaBA) reveals extensive effects of fencing on 
#   wide‐ranging ungulates. J Appl Ecol 58, 690–698 (2021). 
#   https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/1365-2664.13806
#   https://github.com/wx-ecology/BaBA/

# Path straightness is the ratio between the displacement distance and the 
# accumulated step length of a trajectory, ranging from 0 (sinuous) to 1 (straight) 

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

# devtools::install_github("wx-ecology/BaBA", force=TRUE)
# rpath <- "C:/Users/mes473/OneDrive - Cornell University/Documents/R_Packages/BaBA/R"
# source(file.path(rpath, "BaBA_fix.R"))
pacman::p_load(sp, BaBA, here)
source(here("02_scripts", "utilities.R"))
setDataPaths('elephant')
load(procpath('ele.rdata'))

# ******************************************************************************
#                             SET UP DATA TO FIT BABA
# ******************************************************************************

# Data points mutation
points_1h <- ele.df %>% 
  filter(FIXRATE == "1 hour") %>% 
  mutate(Animal.ID = BURST, 
         date=as.POSIXct(DATE.TIME))

points_5h <- ele.df %>% 
  mutate(Animal.ID = BURST, date=as.POSIXct(DATE.TIME)) %>% 
  filter(FIXRATE == '5 hours')

# freq
points_1h %>% nog() %>% 
  mutate(difftime = date - lag(date),
         difftime = ifelse(BURST != lag(BURST), NA, difftime),
         flag = between( as.n(difftime), 0, 200)) %>% 
  group_by(BURST) %>% 
  summarize(n=sum(flag, na.rm=TRUE)) %>% 
  arrange(n)

# getting m/f data and #indivs
mf <- ele.df %>% nog() %>% 
  group_by(ID, SEX) %>% 
  summarize(SEX=first(SEX), NAME = first(ID)) %>% 
  column_to_rownames('ID')
mf$ID <- paste0(mf$SEX, 1:nrow(mf))
nindiv <- ele.df %>% nog() %>% 
  group_by(SEX) %>% 
  summarize(n_distinct(ID))

# ******************************************************************************
#                                 BABA PARAMETERS
# ******************************************************************************
b_time <- 10 #max #hours for behavior = short barrier interaction (bounce or quick)
p_time <- 100 #min #hours for behavior = prolonged barrier interaction (i.e. trapped)
w <- 10 #moving window size, in days.
max_cross <- 3 #max crosses allowed in trace and back-and-forth behavior
d <- 800 #distance considered "fence encounter"

# ******************************************************************************
#                             BARRIER ANALYSIS RUN
# ******************************************************************************
baba_1h <- BaBA(
  animal = points_1h %>% filter(ID==4), 
  barrier = fences, 
  d = d, b_time = b_time, p_time = p_time, w = w, max_cross = max_cross,
  interval = 2, units = "hours", tolerance = 10
)
baba_5h <- BaBA(
  animal = points_5h, barrier = fences, 
  d = d, b_time = b_time, p_time = p_time, w = w, max_cross = max_cross,
  interval = 5, units = "hours", tolerance = 25
)

# combine the two data frames
ambiguous = c("Trace", "Back_n_forth", "Trace_OR_Back_n_Forth")
encounters <- rbind(baba_1h$classification, baba_5h$classification) %>% 
  mutate(SPECIES = mf[AnimalID, "SPECIES"],
         SEX = mf[AnimalID, "SEX"],
         ID = mf[AnimalID, "ID"],
         NAME = mf[AnimalID, "NAME"],
         eventTYPE = ifelse(eventTYPE %in% ambiguous, "Trace", eventTYPE),
         eventTYPE = gsub("_", " ", eventTYPE),
         SEASON = seasons[as.c(date(start_time)),'SEASON'],
         start_side = "",
         end_side = "") %>% 
  relocate(start_side, end_side, .after=end_time) 

# decide what side the beginning and end points are on
# this takes ~2min
nb <- geodata::gadm("Namibia", level=0, path='tmp') %>% st_as_sf()
pb <- txtProgressBar(style = 3)
NE <- nrow(encounters)
for (i in 1:NE) {
  
  # set progress bar and pull matching data
  setTxtProgressBar(pb, i/NE)
  id <- encounters$AnimalID[i]
  data_i <- collar.df %>% filter(COLLAR_ID == id)
  
  # find start and finish points +-1 buffer
  t0 <- encounters$start_time[i]
  t1 <- encounters$end_time[i]
  inx.start <- which(data_i$DATETIME == t0) - 1
  inx.end <- which(data_i$DATETIME == t1) + 1
  
  # then intersection with Namibia to find sides
  ints <- data_i[c(inx.start,inx.end),] %>% 
    st_intersects(nb) %>% as.n()
  sides = ifelse(is.na(ints), "Botswana", "Namibia")
  encounters$start_side[i] <- sides[1]
  encounters$end_side[i] <- sides[2]
}

# assign "jag" and "cross" and rename eventTYPE to EVENT
encounters <- encounters %>% 
  mutate(EVENT = ifelse(eventTYPE == "Cross", 
                        ifelse(start_side == end_side, "Jag", "Cross"), 
                        eventTYPE))

outfile <- here::here("03_output", "barriers", "encounters.cbpp.RData")
save(encounters.cbpp, file=outfile)
# load(outfile)



# ******************************************************************************
#                             CBPP ANALYSIS RUN
# ******************************************************************************
iko <- 
  samo <- cbpp %>% filter(grepl("Samochima", Name)) %>% dplyr::select(geometry)

findEncounters <- function(time, fencename) {
  fencename = toupper(fencename)
  barrier = cbpp %>% filter(grepl(fencename, toupper(Name))) %>% dplyr::select(geometry)
  if (time == "1h") {interval=1; tolerance=10; animal=points_1h
  } else {interval=5; tolerance=25; animal=points_5h;}
  baba <- .BaBA(
    animal = animal, barrier = barrier, 
    d = d, b_time = b_time, p_time = p_time, w = w, max_cross = max_cross,
    interval = interval, units = "hours", tolerance = tolerance
  ) 
  df <- baba$classification %>% as.data.frame() %>% mutate(FENCE=fencename)
  return(df)
}
baba_1h_sam <- findEncounters("1h", "Samochima")
baba_1h_iko <- findEncounters("1h", "Ikoga")
baba_5h_sam <- findEncounters("5h", "Samochima")
baba_5h_iko <- findEncounters("5h", "Ikoga")

# combine the two data frames
ambiguous = c("Trace", "Back_n_forth", "Trace_OR_Back_n_Forth")
encounters.cbpp <- rbind(baba_1h_iko, 
                         baba_5h_sam
) %>% 
  mutate(SPECIES = mf[AnimalID, "SPECIES"],
         SEX = mf[AnimalID, "SEX"],
         ID = mf[AnimalID, "ID"],
         NAME = mf[AnimalID, "NAME"],
         eventTYPE = ifelse(eventTYPE %in% ambiguous, "Trace", eventTYPE),
         eventTYPE = gsub("_", " ", eventTYPE),
         SEASON = seasons[as.c(date(start_time)),'SEASON'])

for (i in 1:nrow(encounters.cbpp)) {
  start <- encounters.cbpp$start_time[i]
  end <- encounters.cbpp$end_time[i]
  df <- collar.df %>%
    filter(COLLAR_ID == encounters.cbpp$AnimalID[i],
           DATETIME %within% interval(start-hours(24), end+hours(24))
    ) %>%
    st_union() %>%
    st_cast("LINESTRING")
  p = ggplot() +
    geom_sf(data=cbpp, color='red') +
    geom_sf(data=df) + 
    ggtitle(i)
  print(p)
}

encounters.cbpp$eventTYPE[c(7,8)] <- "Jag"
encounters.cbpp$eventTYPE[c(13)] <- "Cross"
encounters.cbpp <- encounters.cbpp[-c(1,2,16)]

outfile <- here::here("03_output", "barriers", "encounters_cbpp.RData")
save(encounters.cbpp, file=outfile)
# load(outfile)
