# Vaccination and immigration rates influence raccoon rabies elimination and recolonization in simulated urban-suburban landscapes
### E.M. Beasley & T. Poisot

## **Abstract**

The raccoon variant of the rabies virus (RRV) is managed in the eastern United States and Canada via distribution of oral rabies vaccine (ORV) baits. The goal of ORV distribution is to reach seroprevalence rates (an index of population immunity) of at least 60%, the threshold thought to eliminate RRV. Seroprevalence rates in urban areas rarely reach target levels, predictably leading to rabies outbreaks. However, many urban areas have spent several years rabies-free, aligning with previous work suggesting RRV can be eliminated from urban areas at below-target seroprevalence rates. Using an agent-based model to simulate raccoon populations in urban landscapes, we examined 1) whether RRV can be eliminated at vaccination thresholds below 60%, 2) whether landscapes with below-target vaccination rates are vulnerable to RRV recolonization, and 3) whether the rate and timing of immigration influences elimination and recolonization. Vaccination and immigration rates influenced elimination probability: elimination was more likely and occurred more quickly in landscapes with higher vaccination rates and less likely in landscapes with higher immigration rates. All immigration variables (immigration rate, immigrant disease prevalence, and immigration timing) influenced the probability of recolonization after rabies was eliminated: recolonization was more likely in landscapes with high immigration rates and when immigrants had higher disease prevalence, but less likely when immigration occurred seasonally rather than continuously. Vaccination did not have a clear effect on recolonization probability but reduced the number of rabies cases during a recolonization event. Although elimination was highly likely in our simulated landscapes due to their small spatial extent, our results suggest that vaccination rates of at least 50% result in timely rabies elimination (median 1.5 years). After elimination is achieved, strategies for preventing infected individuals from entering the rabies-free area are essential for preventing recolonization events, as vaccination rates had a much smaller effect on recolonization that immigration rates and timing. Understanding long-distance movements of host individuals is crucial for managing diseases such as rabies which likely persist at the regional scale.

## **Data**

### BurlingtonLandCover2016.grd, BurlingtonLandCover2016.gri

NLCD land cover data from the greater Burlington, Vermont area at a 30m x 30m resolution, clipped to an extent of -73.347 – -73.0.17 degrees longitude and 44.374 – 44.587 degrees latitude.

## **Code**

Scripts for the main analysis are intended to be executed in the following order.

### LandCoverPropCalculations.R

Calculates land cover proportions and spatial autocorrelation of land cover data in the greater Burlington, VT area.

### Functions.jl

Defines functions needed for the agent-based model described in the manuscript.

### CreateParams.jl

Creates a .csv file containing all possible parameter combinations described in the manuscript, and writes a .sh file for executing the array job on the Narval computing cluster hosted by the Digital Research Alliance of Canada.

### RunTheSim.jl

Executes the agent-based model described in the manuscript and saves the results as a .csv file.

### figures_full.R

Processes agent-based model outputs and visualizes the results.

## **Subfolders**

### outputs

Contains all results of the agent-based model described in the main manuscript saved in .csv format.

### FunctionalityTest

Contains files for testing various functions of the agent-based model, specifically landscape construction, agent movement, transitions between disease states, and agent mortality. A README file with more details can be found in the subfolder.

### ParamSensitivity (CURRENTLY UPDATING)