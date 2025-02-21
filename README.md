# challenging_neuromarkers
Here, we provide the corresponding R code to the simulation our paper 'The troubling pursuit of brain biomarkers for mental health: multiple realizability, emergence and noise.' is based on. With the code, our findings can be replicated, or the reader may experiment with different parameters to see how these impact diagnostic model parameters.

The code consists of two parts:

Part 1 - Technical Script ('functions'): In this script, we load necessary packages, then pre-define the simulation and plot functions, that are called in the main script. Mathematical and technical definitions can be found in this part of the code.

Part 2 - Practical Script ('simulation'): In this script, we define the variables of interest and simulate the data. We call functions that were defined in the technical script. This script allows the user to test different model parameters (e.g., different vectors of brain-symptom correlations) and see how changes in parameters impact diagnostic model parameters such as R squared or classification accuracy, both, on an individual or group level. 
