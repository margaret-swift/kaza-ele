# KAZA Elephant Behavior Modeling Project
AHEAD x WWF collaboration with funding from Atkinson Center for Sustainability

Margaret Swift, Robin Naidoo, Steve Osofsky, Shirley Atkinson, Anna Songhurst

**WARNING: This repository is public, but raw data files (especially elephant GPS data) are sensitive and never stored or tracked through GitHub.**

Agent-based model code is housed at: [abmAME](https://github.com/margaret-swift/abmAME) 

| Step | Progress |
| -----| -------- |
| [Gather Data](#1gather-data)| 100% |
| [EDA](#2-exploratory-data-analysis) | 100% |
| [Simulate basic movements](#3simulating-elephant-movements) | 100% |
| [Model landscape use](#4modeling-landscape-use) | 0% |
| [Validation and improvement](#5validation) | 0% |
| [Display](#6displaying-simulations) | 0% |
| [Extensions](#7extending-simulations) | 0% |

### Tables
- [Spatial data](#spatial-data)
- [Patterns and metrics](#patterns-and-metrics)
- [Elephant step statistics](#elephant-step-statistics)

# Project steps
## 1.	Gather data
  - [x] **Data for elephants** collared in Namibia have been provided by **Robin Naidoo**; I have reached out to **Anna Songhurst** (11/6/23) about a possible collaboration to bring her expertise and data for elephants collared in Botswana.
  - [x] **Landscape-level spatial data** have been collected and are provided [below](#spatial-data).
        
_[^Top^](#kaza-elephant-behavior-modeling-project)_

## 2. Exploratory Data Analysis
  - [x] **Agent movement statistics** have been collected in [02_scripts/2_eda/eleStats.R](https://github.com/margaret-swift/kaza-ele/blob/main/02_scripts/2_eda/eleStats.R); these data are provided [below](#elephant-step-statistics), and include:
    - [x] Home range size (female only), by **season**
    - [x] Rates of fence crossing by **activity type** and **season** (male only)
    - [x] Step size by **state, sex, season**
  - [x] **Activity states** have been defined using a Hidden Markov Model ([McClintock & Michelot 2018](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12995)). Also approached was the M4 Model ([Cullen et al 2021](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13745)); I may return to this later.
    - [x] I need to re-run this model with **season** as well as **sex** and gather transition matrices.
        
_[^Top^](#kaza-elephant-behavior-modeling-project)_

## 3.	Simulating basic elephant movements 
  - [x] **Basic agent-based model** has been implemented. At first, I tried [SimRiv](https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-019-0154-8), but found that its structure did not allow for flexibility in seasonal/diel cycles or attractive landmarks (waterholes). I then chose [abmAnimalMovement](https://f1000research.com/articles/11-1182) and found that this structure was flexible enough for our needs.
  - [x] **Augment chosen ABM** Currently, I am in the process of augmenting the abmAnimalMovement code in a new package I've termed '[abmFences](https://github.com/margaret-swift/abmFences)'. The following is a list of important features that need to be added to the code in order to move forward:
    - [x] ~~Fence response behavior~~
    - [x] ~~Selectively permeable fences with probability changing based on sex, season, and activity state~~
    - [x] ~~Attraction to range centroids~~
    - [x] ~~Interactions with landscape raster values~~

_[^Top^](#kaza-elephant-behavior-modeling-project)_

## 4.	Modeling Landscape Use 
  - [ ] **Integrated Step Selection**: In order to more accurately simulate elephant movements, we have to understand how they use the landscape currently and then transfer this knowledge onto the spatial data to calculate accurate resistance rasters. To do this, we will run an Integrated Step Selection Function ([iSSF](https://www.biorxiv.org/content/10.1101/2023.08.10.552754v1)) to estimate elephant responses to various landscape features, then apply the results of this model in the next step.
  - [ ] **Create resistance rasters**: Spatial data should be transformed into rasters that represent an agent's willingness to travel through each cell, depending on **sex, season,** and **activity state.** These rasters should hold values from 0 to 1, where higher values are more likely to be chosen (see [Marshall and Duthie 2022](https://f1000research.com/articles/11-1182), Fig. 3).
        
_[^Top^](#kaza-elephant-behavior-modeling-project)_

## 5. Validation
  - [ ] **Define metrics**: Define and link specific quantitative metrics (see [Butts et al 2022](https://www.sciencedirect.com/science/article/pii/S0304380022001132)) to the qualitative patterns we think are necessary to replicate (table below), using EDA to explicitly define these characteristics. Here we might also define different movement or activity states, depending on the modeling method. It might also be a good idea to run a Barrier Behavior Analysis (BaBA, Xu et al 2021) on the elephant data we choose to use, so we can then run the same analysis on the simulated data & see if the encounter behavior is similar.
  - [ ] **Model selection**: Calculate metrics from step 1 and compare models. Which ones accurately replicate the critical spatial patterns? Is any homebrewing necessary?
  - [ ] **Model improvement**: Temporal and individual elements should then be added to the model to better reflect the seasonality of water supply, vegetation, temperature, and elephant movements. In addition, individual elephants likely would remember waterholes and fence gaps in particular; they shouldn’t be treated like a random molecule. In addition, bulls are known to create fence gaps; should agents be able to modify their landscape in this capacity?

_[^Top^](#kaza-elephant-behavior-modeling-project)_

## 6.	Displaying simulations
  - [ ] **Model with and without key fences**
  - [ ] **Model with different climate scenarios**
  - [ ] **Create display of our results**
        
_[^Top^](#kaza-elephant-behavior-modeling-project)_

## 7.	Extending simulations
  - [ ] **Activity states**: Could draw from TOD, time since last state change, etc. and determine whether the agent is foraging, exploring, or resting. This would then determine the characteristics of the distribution from which we’re drawing step lengths and turning angles (internal state x landscape resistance perhaps). 
  - [ ] **Additional Species**: Once we have the model working, we could extend this to any number of species, including roan, oryx, buffalo, or even a general species X with Y characteristics. Steve suggests prioritizing buffalo and cattle since they’re the impetus for putting these fences up in the first place. Robin suggests using buffalo data from the Caprivi strip, Shirley on zebra; zebra and buffalo don’t cross often. What about expanding more beyond N and B into the rest of KAZA? 
  - [ ] **Fence Structure**: Right now we’ll just treat all the fences in our dataset as the same type, but this isn’t quite the case. Is there a way to verify which fences are robust versus falling apart? Or perhaps if an extreme fence were replaced with a smaller barrier that cattle still couldn’t cross (i.e. lower the cost of crossing for elephants).
  - [ ] **Future Infrastructure**: E.g., the railway across the Caprivi strip (along the highway) or if they paved certain corridors. Railway in NE Hwange. Shirley mentioned that some of the waterholes in Hwange aren’t being filled due to budget constraints, and elephants are moving into Botswana earlier in the year to get water. What might the removal of key waterholes in Hwange or Khaudum do to elephant movements (and fence encounters)? 

_[^Top^](#kaza-elephant-behavior-modeling-project)_

# Tables
## Spatial data
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

_[^Top^](#kaza-elephant-behavior-modeling-project)_

## Patterns and Metrics
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

_[^Top^](#kaza-elephant-behavior-modeling-project)_

## Elephant step statistics
See [Butts et al 2022](https://www.sciencedirect.com/science/article/pii/S0304380022001132) for a more detailed look at this EDA approach for creating an ABM.

### Feature permeability (literature)
From [Naidoo et al 2022](https://www.frontiersin.org/articles/10.3389/fcosc.2022.788133/).
Percentage of crossings @ 1km encounter threshold
(should we use 25km instead? What is the utility of a larger threshold if we are controlling the movements?)

| | River	| Road	| Fence| 
| --- | -------| --------| ------| 
| Female	| 10.1	| 15.3 | 0| 
| Male	| 14.5	| 25.8	| 3.5| 

_[^Top^](#kaza-elephant-behavior-modeling-project)_

### Average step size (m) and speed (meters per second)
| Sex | Activity state | Season | distance $\mu$ | distance $\sigma$ | speed $\mu$ | speed $\sigma$ |
| --- | -------------- | ------ | -------------- | ----------------- | ----------- | -------------- |
| F | correlated walk |DRY  |  671 |    829 |0.183    |0.240|
| F | correlated walk |WET  |  777 |    883 | 0.206   |0.240|
| M | correlated walk | DRY |  598 |    854 | 0.152   |0.220|
| M | correlated walk |WET  |  559 |    692 | 0.154   |0.197|
| | | | | | | |
| F | foraging |    DRY |  563|    734| 0.164     |0.257|
| F | foraging |    WET |  630|    817| 0.170     |0.217|
| M | foraging |    DRY |  554|    791| 0.138     |0.187|
| M | foraging |    WET |  451|    581| 0.124     |0.157|
| | | | | | | |
| F | resting |     DRY |  524|    737| 0.144     |0.225|
| F | resting |     WET |  620|    838| 0.164     |0.221|
| M | resting |     DRY |  486|    717| 0.124     |0.173|
| M | resting |     WET |  468|    581| 0.131     |0.182|

_[^Top^](#kaza-elephant-behavior-modeling-project)_

### Average home range size (females only)
| Season | $\mu$ (m2) | $\sigma$ |
|--------| ----- | ---------|
| Dry	season | 191,738 | 260,010 |
| Wet season | 282,596 | 250,076 |

_[^Top^](#kaza-elephant-behavior-modeling-project)_

### Average barrier crossing rates

#### Fences (males only)
| Season | Activity state | Value | 
|--------| -------------- | ---- |
| Dry season | Resting | XX |
| Wet	season | Resting | XX |
| Dry season | Foraging | XX |
| Wet	season | Foraging | XX |
| Dry season | Exploring | XX |
| Wet season | Exploring | XX |

#### Rivers
| Season | Sex | Value | 
|--------| --- | ----- |
| Dry season | Female | XX |
| Dry season | Male | XX |
| Wet	season | Female | XX |
| Wet	season | Male | XX |

#### Roads
| Season | Sex | Value | 
|--------| --- | ----- |
| Dry season | Female | XX |
| Dry season | Male | XX |
| Wet	season | Female | XX |
| Wet	season | Male | XX |

_[^Top^](#kaza-elephant-behavior-modeling-project)_
