# Character-Displacement-Storage-Effect
Here is a collection of code to run a evolutionary model of germination patterns in two species of annual plants in a fluctutaing environment. The species coexist via the temporal storage effect.

There are 2 functions that are required for the analysis. 

["Ecological_Dynamics.m"]
A function that takes in model parameters and outputs total community stabilizing mechanisms (A\bar) and adjusted species average fitness difference (\kappa_i\prime and \kappa_j\prime) exactly from invader growth rates and also gives approximations under an assumption of small varaibility in environmental parameters. 

Inputs:
  1. thetag = a 1x2 vector of \theta){G_i} values for each species i
  2. thetav = a 1x2 vector of \theta_{V_i} values for each species i
  3. gen = the number of time steps to evaluate the long-term growth rate
  4. bunrin = the number of time steps to omit in order to minimize effects of transient dynamics on the average long-term growth rate
  5. Parameters = a cell array with 8 elements
    a. s = seed survival fraction
    b. y = a 1x2 vector of species low density yield values
    c. rho = the information content of the cue
    d. muG = mean EG
    e. muV = mean EV
    f. sigmaG = standard deviation of EG
    g. sigmaV = standard deviation of EV
    h. alpha = the competition coefficient
    
Outputs:
  1. A = (r_1 + r_2)/2; exact (within numerical error) community average stabilizing mechanisms from simulation
  2. SpecAveFit = (r_1 - r_2)/2; adjusted species average fitness difference (\kappa_1\prime - \kappa_2\prime) from low density growth rates computed by simulation.
  3. kappaDiff = \kappa_1 - \kappa_2 from the small varaince approximation. 
  4. ComAveDeltaJ = community average relative nonlinearity from small variance approximation.
  5. ComAveDeltaIG = community average storage effect from germination, calculated from the small variance approximation.
  6. ComAveDeltaIV = community average storage effect from vigor, calculated using the small variance approximation. 

Dependencies: None


["AdaptDynamics22_PolarCoord.m"]
A function that runs the adaptive dynamics simulations under evolution of both competitors. The output here is are the dynamics of the species traits under the assumptions of adaptive dynamics.

Inputs: 
  1. thetaginit = the initial germination trait values, (\theta_{G_1}, \theta_{G_2))
  2. thetav = the values of vigor traits for each species (\theta_{V_1}, \theta_{V_2})
  3. gen = number of timesteps to evaluate mutant long-term growth rate
  4. mutnum = number of mutations to run. The number of mutations is split evenly between the two species such that species 1 has mutnum/2 mutations during the simulation.
  5. Parameters = a cell array with 8 elements
    a. s = seed survival fraction
    b. y = a 1x2 vector of species low density yield values
    c. rho = the information content of the cue
    d. muG = mean EG
    e. muV = mean EV
    f. sigmaG = standard deviation of EG
    g. sigmaV = standard deviation of EV
    h. alpha = the competition coefficient

Outputs:
  1. a 1xmutnum vector of \theta_{G_1} values for each mutational time step
  2. a 1xmutnum vector of \theta_{G_2} values for each mutational time step
  3. a 2xmutnum matrix of values of E(N_1) and E(N_2) for each mutational time step. This output can be used to see if the species population reaches extremely low population sizes during its evolutionary timecourse. 

Dependencies: None

The other files in this repository are:
["viridis.m"] 
A colormap used in plotting.  ** This is not my file. It was acquired open access from https://bids.github.io/colormap/ **

["Selection_Portraits.m]
A script to generate the selection portrait figures used to understand long-term evolutionary equilibria. 
This script has three examples: the first is to visualize the effect of changes in seed yield, the second is to visualize the effect of competitive asymmetries, and the third is to visualize the effect of species differences in vigor. 
The script is broken into these three sections and each can be run without running the others by specifying a 3-element vector at the beginning of the script. The vector is called "RUN" and each element takes values 0 or 1. If the element is 1, then that section of the script is run. If it is not 1, the section is not run. 

Dependencies: "AdaptDynamics22_PolarCoord.m"; "Ecological_Dynamics.m"; "viridis.m"

["ESC_FULL.m]
A script to find evolutionary endpoints across a range of parameters and the associated coexistence quantities Abar and adjusted species average fitness differences. 
The script is broken into four sections: a section changing the value y and sigmaV, all esle fixed; a section changing the difference in vigor traits for different values of rho, all esle fixed; a section changing the difference lny1 - lny2 for different values of rho, all esle fixed; a section plotting the effect of rho for different values of sigma_{E_V}, all else fixed.
Each section can be run without running the others by specifying a 4-element vector "RUN" at the beginning of the script. Setting any element of this vector to 1 runs that section. If the value is not 1, the section is not run. 

Dependencies: "AdaptDynamics22_PolarCoord.m"; "Ecological_Dynamics.m"; "viridis.m"

