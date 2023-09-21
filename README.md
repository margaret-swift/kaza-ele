# KAZA Elephant Behavior Modeling Project
AHEAD x WWF collaboration with funding from Atkinson Center for Sustainability

Margaret Swift, Robin Naidoo, Steve Osofsky, Shirley Atkinson, Martin Gilbert, Helen Lee

# Project steps
## 1.	**Agent data**
  - [ ] **Gather elephant data**: Robin should be sending elephant data over soon.
  - [ ] **Gather agent attributes**: List movement characteristics of agents from data and gather missing parameters from literature. This includes a deeper dive into literature on elephant preferences for certain landcover types (although we can find that out through modeling our data) or responses to settlement areas, for example.
 
## 2.  **Spatial data**
  - [x] ~~**Gather spatial data**: List spatial features and gather spatial data~~. DONE 9/21/23
  - [ ] **Combine spatial data**: Depending on the model, there should be two separate layers: One “resistance” layer representing how difficult the landscape is to move through (landcover types, rasterized linear features, maybe Robin’s Circuitscape output?), and one “magnet” layer with features that either attract or repel elephants (water features, urban areas, cropland).
  
## 3.  **Exploratory Data Analysis**
  - [ ] **Define metrics**: Define and link specific quantitative metrics (see Butts et al 2020) to the qualitative patterns we think are necessary to replicate (table below), using EDA to explicitly define these characteristics. Here we might also define different movement or activity states, depending on the modeling method. It might also be a good idea to run a Barrier Behavior Analysis (BaBA, Xu et al 2021) on the elephant data we choose to use, so we can then run the same analysis on the simulated data & see if the encounter behavior is similar.

## 4.	**Modeling**
  - [ ] **Basic model implementation**: Run a basic version of SimRiv and abmMovement with agent parameters from (1b) and simplified spatial rasters. This way we can make sure the behaviors seem right before dedicating the computational power needed for the main, large raster.
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
 
