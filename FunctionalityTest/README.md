# /FunctionalityTest
### Vaccination and immigration rates influence raccoon rabies elimination and recolonization in simulated urban-suburban landscapes
### E.M. Beasley & T. Poisot

Contains code for testing various functions of the agent-based model, specifically landscape construction, agent movement, transitions between disease states, and agent mortality.

## **Code**

### Functions_smol.jl

Defines functions needed to test the agent-based model. Many functions are modified from those in the main text for use with a smaller simulated landscape, which reduces computation time. Other changes from the functions in the main text create and save simulated data for determining whether the function is operating as expected.

### FunctionalityTests.jl

Executes functionality tests of the agent-based model and saves outputs as .csv files. The user can define which component to test by commenting/uncommenting certain sections of code.

### FunctionalityTests.R

R script for generating figures of the results of the functionality tests.