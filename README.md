## Computer Code for "Mechanisms underlying higher order interactions: from definitions to ecological processes"

This directory contains all R scripts necessary to run the simulations and analyses and generate the figures presented in the manuscript. 

### Reproducing analyses 

The 'code' directory contains all the R scripts. To recreate the analyses in the manuscript open an R session and set the working directory to this folder containing this README file (alternatively open the .Rproject file in Rstudio).  Then run the "run_all_scripts.R" script. **Note: running all analyses may take 15 to 20 minutes to complete.**

### File details 

1. run_all_scripts.R

  Make a figures and output folder in the working directory and run all the scripts in order above. 

2. code/set_up_parms.R 

  Defines the simulation parameters used by the mechanistic growth model (Table S1). Generates figure 4. 

3. code/run_simulations.R

  Iterate through all the combinations of competitor density explored in the main text. 
  
4. code/fit_models.R 

  Fit the phenomenological competition models to the simulated experimental data generated by the run_simulations.R script. 

5. code/plot_fits.R
  
  Generate figure 5, and 6 in the main text and figure S2 in the supporting information. 
  
6. code/plot_parameters.R

  Generate figure 7 in the main text. 

7. code/vary_tradeoff.R 

  Complete simulation for Appendix A "The effect of trait differences on higher order interactions". Vary the parameters defining the three species and then compare the changes in higher order interaction strength. Takes several minutes to run.  
  
8. code/plot_trade_off_results.R

  Generate figure comparing species trait differences to the strength of higher order interactions. Generate Figures A1 and A2 in Appendix A.  

9. code/simulation_functions.R 

  The differential equation functions used by the simulation. 

10. code/process_sim_data.R 

  Function to organize the results of the simulations. 
  
11. code/phenomenological_models.R 
  
  Functions for each phenomenological model used to fit the simulation data. 
  
12. code/plotting_parameters.R
  
  Define variables (colors etc.) used to generate figures. 
  
### Built With 

platform       x86_64-apple-darwin15.6.0   
arch           x86_64                      
os             darwin15.6.0                
system         x86_64, darwin15.6.0        
status                                     
major          3                           
minor          6.1                         
year           2019                        
month          07                          
day            05                          
svn rev        76782                       
language       R                           
version.string R version 3.6.1 (2019-07-05)
nickname       Action of the Toes   

### Required R packages 

1. tidyverse_1.2.1
2. gridExtra_2.3   
3. deSolve_1.24   
4. stringr_1.4.0  

### Authors 

Obscured for peer review 






