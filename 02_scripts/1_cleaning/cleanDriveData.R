# cleanDriveData
# Created 6 Feb 2024
# Margaret Swift <margaret.swift@cornell.edu>

# ******************************************************************************
#                             DATA & LIBRARY LOADING
# ******************************************************************************

here::i_am('02_scripts/1_cleaning/cleanDriveData.R')
source(here::here('02_scripts', 'utilities.R'))
pacman::p_load(googledrive, patchwork, nlme)
setDataPaths('elephant')
load(procpath('ele.rdata'))

drivepath = "GEEData/ele_covariates"
files = drive_ls(drivepath, pattern=".csv$")

drive.df <- NULL
for (i in 1:nrow(files)) {
  fname <- files[i,]
  message(gsub('.csv', '', fname$name))
  f <- drive_read_string(fname, encoding="UTF-8") %>% read.csv(text=.)
  if (is.null(drive.df)) drive.df <- f
  else drive.df <- rbind(drive.df, f)
}

ele.khau <- left_join(ele.df, drive.df, by="INX") %>% 
  left_join( pdata.f, by="INX") %>% 
  st_intersection(khau) %>% 
  nog() %>% 
  filter(!is.na(STATE), year(DATE) > 2013)
names(ele.khau) <- gsub('.x', '', names(ele.khau))
inx.rm <- grepl('\\.y', names(ele.khau))
ele.khau <- ele.khau[,!inx.rm]


save(ele, ele.df, ele.khau, drive.df, file=procpath('eleKhau.rdata'))
