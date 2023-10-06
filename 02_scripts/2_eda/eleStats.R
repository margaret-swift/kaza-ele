# eleStats.R
# Created 04 Oct 2023
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/2_eda/eleStats.R')
source(here::here('02_scripts','utilities.R'))
pacman::p_load(sp, adehabitatHR, reshape2)
setDataPaths('elephant')
load(here(procpath, 'ele.rdata'))


# ******************************************************************************
#                           General path statistics
# ******************************************************************************

# sort individuals by # days collared
ele.df %>% nog() %>% 
  group_by(ID) %>% 
  summarize(min=min(DATE),
            max=max(DATE),
            days = max-min) %>% 
  arrange(days) %>% 
  print(n=43)

# sort individuals by fixrate type and give count
ele.df %>% nog() %>% 
  group_by(ID, FIXRATE) %>% 
  summarize(n=n()) %>% 
  arrange(FIXRATE) %>% 
  print(n=120)

# average speeds
sums <- ele.df %>% nog() %>% 
  group_by(SEX) %>% 
  filter(FIXRATE == "1 hour") %>% 
  summarize(MPS.m = mean(MPS, na.rm=TRUE),
            MPS.med = median(MPS, na.rm=TRUE),
            MPS.sd = sd(MPS, na.rm=TRUE),
            MPS.max = max(MPS, na.rm=TRUE),
            DIST.m = mean(DIST, na.rm=TRUE),
            DIST.med = median(DIST, na.rm=TRUE),
            DIST.sd = sd(DIST, na.rm=TRUE),
            DIST.max = max(DIST, na.rm=TRUE),
            n=n())

d <- ele.df %>% nog() %>% 
  filter(FIXRATE=="1 hour") %>% 
  mutate(SEX=factor(SEX, levels=c("M", "F")))
ggplot(d, mapping=aes(x=DIST, fill=SEX)) + 
  geom_density(position="dodge", bins=100, alpha=0.4) + 
  geom_vline(data=sums, 
             mapping=aes(xintercept=DIST.m, color=SEX)) + 
  xlim(0, 3000) + 
  scale_fill_brewer(palette="Accent") + 
  scale_color_brewer(palette="Accent") + 
  plot.theme

# ******************************************************************************
#                           MATCHING PATTERNS TO STATISTICS
# ******************************************************************************

# Elephant movements trace along fences and channel along omiramba



# Pinch points across roads and at rivers	



# Habitual/repeated movement to and from water sources, especially artificial waterholes		



# Deflection/permeability differences by sex and boundary type (fence [and fence type], road, river)



# Elephants move away from water sources in the wet season (expansion contraction)



# 80% of roan crossings were during the wet season. Does this happen for elephant too?	



# Elephant attraction to some areas with higher quality resources?	