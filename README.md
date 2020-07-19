# PogZ_paper
 
## Analysis of LFP experiments, including LFP data, annotating relevant behavior data, visualization

- pairElectrodeBlockAnalysis is master function for calculating power & lag lead statistics between 2 electrodes averaged over blocks of time (eg time in open arm vs time in closed arm). pairElectrodeTransitionAnalysis is a mster function for calculating power & lag lead statistics between 2 electrodes around a single time point (eg point at which a mouse enters an open arm)

- pairElectrodeBlockAnalysis and pairElectrodeTransitionAnalysis call expPower and expLagLead for individual calculations for each file. expLagLeag calls trialLagLead for individual block measurements

- annotationFunctions.py takes position data from anymaze output and converts it into annotation files for the LFP analysis. Block analysis requires 3 columns - start time, end time, and run type. Transition analysis is two columns  - timepoint of interest and run type.

## Analysis of patch clamp experiments
- allCellProperties is a master funciton for calculating intrinsic cell properties from CCIV recordings
- optoESPC, optoIPSC, and optoSpikes are functions for calculating the response to light flashes in optogenetic experiments.
