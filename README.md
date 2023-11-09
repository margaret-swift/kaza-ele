# KAZA Elephant Behavior Modeling Project
AHEAD x WWF collaboration with funding from Atkinson Center for Sustainability

Margaret Swift, Robin Naidoo, Steve Osofsky, Shirley Atkinson, ...

**WARNING: This repository is public, but raw data files (especially elephant GPS data) are sensitive and never stored or tracked through GitHub.**

Agent-based model code is housed at: www.github.com/margaret-swift/abmFences


# Project steps
## 1.	**Agent data**
  - [x] ~~**Gather elephant data**: Robin should be sending elephant data over soon.~~. DONE 9/16/23
  - [x] ~~**Gather agent attributes**: List movement characteristics of agents from data and gather missing parameters from literature. This includes a deeper dive into literature on elephant preferences for certain landcover types (although we can find that out through modeling our data) or responses to settlement areas, for example.~~ DONE, saved in 02_scripts/2_eda/eleStats.R
    - [ ] Definition and transition matrix for activity states by **sex** and **season**
    - [ ] Home range size (female only), by **season**
    - [ ] Rates of fence crossing by **activity type** and **season** (male only)
    - [ ] Step size by **state, sex, season**
 
## 2.  **Spatial data**
  - [x] ~~**Gather spatial data**: List spatial features and gather spatial data~~. DONE 9/21/23
  - [x] **Combine spatial data**: Depending on the model, there should be two separate layers: One “resistance” layer representing how difficult the landscape is to move through (landcover types, rasterized linear features, maybe Robin’s Circuitscape output?), and one “magnet” layer with features that either attract or repel elephants (water features, urban areas, cropland).
  
## 3.  **Exploratory Data Analysis**
  - [ ] **Define metrics**: Define and link specific quantitative metrics (see [Butts et al 2022](https://www.sciencedirect.com/science/article/pii/S0304380022001132)) to the qualitative patterns we think are necessary to replicate (table below), using EDA to explicitly define these characteristics. Here we might also define different movement or activity states, depending on the modeling method. It might also be a good idea to run a Barrier Behavior Analysis (BaBA, Xu et al 2021) on the elephant data we choose to use, so we can then run the same analysis on the simulated data & see if the encounter behavior is similar.
  - [x] ~~**Segment activity states**: Segment GPS paths by activity state~~
    - [x] ~~Hidden Markov Model ([McClintock & Michelot 2018](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12995))~~
    - [x] ~~M4 Model ([Cullen et al 2021](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13745))~~

## 4.	**Modeling**
  - [ ] **Basic model implementation**: Run a basic version of [SimRiv](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-019-0154-8) and [abmMovement](https://f1000research.com/articles/11-1182) with agent parameters from (1b) and simplified spatial rasters. This way we can make sure the behaviors seem right before dedicating the computational power needed for the main, large raster.
    - [x] ~~Fence response behavior~~
    - [x] ~~Selectively permeable fences with probability changing based on sex, season, and activity state~~
    - [x] ~~Attraction to range centroids~~
    - [ ] interactions with landscape raster values
  - [ ] **Basic model assessment**: With controllable, small-scale environmental rasters, how do the movements we are seeing in the agents matching up against the metrics determined in 3a?
  - [ ] **Realistic model implementation**: Once the basic model is up and running, add in spatial layers from (2b).
  - [ ] **Model selection**: Calculate metrics from step (3a) and compare models. Which ones accurately replicate the critical spatial patterns? Is any homebrewing necessary?
  - [ ] **Model improvement**: Temporal and individual elements should then be added to the model to better reflect the seasonality of water supply, vegetation, temperature, and elephant movements. In addition, individual elephants likely would remember waterholes and fence gaps in particular; they shouldn’t be treated like a random molecule. In addition, bulls are known to create fence gaps; should agents be able to modify their landscape in this capacity?

## 5.	**Display**
  - [ ] **Model with and without key fences**
  - [ ] **Create display of our results**

## 6.	**Extensions**
  - [ ] **Activity states**: Could draw from TOD, time since last state change, etc. and determine whether the agent is foraging, exploring, or resting. This would then determine the characteristics of the distribution from which we’re drawing step lengths and turning angles (internal state x landscape resistance perhaps). 
  - [ ] **Additional Species**: Once we have the model working, we could extend this to any number of species, including roan, oryx, buffalo, or even a general species X with Y characteristics. Steve suggests prioritizing buffalo and cattle since they’re the impetus for putting these fences up in the first place. Robin suggests using buffalo data from the Caprivi strip, Shirley on zebra; zebra and buffalo don’t cross often. What about expanding more beyond N and B into the rest of KAZA? 
  - [ ] **Fence Structure**: Right now we’ll just treat all the fences in our dataset as the same type, but this isn’t quite the case. Is there a way to verify which fences are robust versus falling apart? Or perhaps if an extreme fence were replaced with a smaller barrier that cattle still couldn’t cross (i.e. lower the cost of crossing for elephants).
  - [ ] **Future Infrastructure**: E.g., the railway across the Caprivi strip (along the highway) or if they paved certain corridors. Railway in NE Hwange. Shirley mentioned that some of the waterholes in Hwange aren’t being filled due to budget constraints, and elephants are moving into Botswana earlier in the year to get water. What might the removal of key waterholes in Hwange or Khaudum do to elephant movements (and fence encounters)? 
 

# Patterns and Metrics
See [Butts et al 2022](https://www.sciencedirect.com/science/article/pii/S0304380022001132) for a more detailed look at this EDA approach for creating an ABM.

| Pattern | Scale | Features  | Metric |
| ------------- |-------------| -----| -------|
| Elephant movements trace along fences and channel along omiramba	| Spatial	| Fences & omiramba | |
| Pinch points across roads and at rivers	| Spatial	| Roads and rivers	| | 
| Habitual/repeated movement to and from water sources, especially artificial waterholes	| Spatial and Temporal	| Natural & artificial waterholes	| | 
| Deflection/permeability differences by sex and boundary type (fence [and fence type], road, river)	| Spatial	| Fence, road, river	| Encounter and crossing rates should be similar to those found in Naidoo et al 2022 (see below) | 
| Elephants move away from water sources in the wet season (expansion contraction)| 	Temporal| 	fences	| | 
| 80% of roan crossings were during the wet season. Does this happen for elephant too? | Temporal | fences ||
| Elephant attraction to some areas with higher quality resources?	| Spatial and Temporal?	| Landcover type?	| | 


# Spatial data

| Description	| Source	| Res| Extent	| Type	| Status| 
| ------------|---------|----|--------|-------|-------|
| Elevation (DEM) and slope (calculated) | [USGS SRTM](https://earthexplorer.usgs.gov/) | 30m	| KAZA 	| Raster	| Obtained, reprojected, and mosaicked on ArcGIS; slope calculated in ArcGIS| 
| Human settlements	| [ESA World Settlement Footprint 2019](https://geoservice.dlr.de/web/maps/eoc:wsf2019) (Sentinel-1 and -2)  | 10m	| Africa	| Raster	| Obtained | 
| Landcover	| WWF (link? What dataset did these come from?)	| 10m	| Africa	| Raster	| Obtained| 
| Ephemeral surface water	| [Schaffer-Smith et al 2022](https://doi.org/10.4211/hs.6f5b34803dc247e890925d7f26b04a3b) (Sentinel-2)  | 10m	| Khaudum 	and Bwabwata| Raster	| Obtained | 
| Fire Incidence | MODIS Burned Area product ([Fire CCI](https://climate.esa.int/en/projects/fire/)) | 250m | KAZA | Raster | Obtained | 
| | | | | | | 					
| Fences	  | Robin... where did he get these?	| NA	| KAZA	| Vector	| Obtained| 
| Roads	    | OpenStreetMap (Angela) (see metadata file)	| NA	| KAZA 	| Vector	| Obtained| 
| Rivers	  | Digitized by Robin from various sources ([GAIA](http://gaia.geosci.unc.edu/rivers/)) | NA	| KAZA	| Vector	| Obtained| 


# Feature permeability
From [Naidoo et al 2022](https://www.frontiersin.org/articles/10.3389/fcosc.2022.788133/).
Percentage of crossings @ 1km encounter threshold
(should we use 25km instead? What is the utility of a larger threshold if we are controlling the movements?)

| | River	| Road	| Fence| 
| --- | -------| --------| ------| 
| Female	| 10.1	| 15.3 | 0| 
| Male	| 14.5	| 25.8	| 3.5| 
