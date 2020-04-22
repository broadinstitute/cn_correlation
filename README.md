# cn_correlation

## History
This algorithm was created by Travis Zack under the direction of Rameen Beroukhim. It was applied to a 5000 patient sampe pan-cancer cohort and published in [Zack213] (Zack et al., *Nat Gen*, 2013). Travis's working logs were adapted post-publication into a re-useable functional format by Steve Schumacher. The "lockstep tempering" variant was invented to improve the effectiveness of the algorithm with a smaller (700ish), single-disease cohort for [Schumacher2017] (Schumacher et al., *PlosOne*, 2017). The original implementation was adapted to use the GridEngine MPE (Broad UGER) and use less memory (allowing more permutations) for a 11K patient sample analysis published in [Gong2016] (Gong et al., *Neoplasia*, 2017). The algorithm has adapted for studies of correlation of arm-level changes in low grade gliomas (TCGA et al., *NEJM*, 2015) and aneuploidy (Taylor et al., *Cancer Cell*, 2018).

This code was originally maintained in the Broad Institute "CancerGenomeAnalysis" Subversion 1.6 repository under the path
https://svnrepos/CancerGenomeAnalysis/trunk/matlab/snp/correlation. It generally depends on some of the code under https://svnrepos/CancerGenomeAnalysis/trunk/matlab/

## Overview
-This method creates a background model for accessing the significance of observed co-occurring (or mutually excluding) copy number events by creating a large number of permutations from the observed data. 

## canonical reference files
- *Binary_amps.mat*
- *Binary_dels.mat*
- *new_samples.mat*
- *D.mat*
- *peak_regs.mat*




## Code File Descriptions
-README.txt - old documentation, more thinking out how to functionally decompose Travis's scripts, less useful for final code.

### Permutation Functions

#### high level script examples
- *ng_pancan_corrperm_hl_recap.m* - script to recapitulate [Zack2013] "high-level" analysis using refactored functions
- *ng_pancan_corrperm_ll_recap.m* - script to recapitulate [Zack2013] "low-level" analysis using refactored functions

#### simulated annealing algorithm
- *corrperm_cooldown_ratio.m* - run one iteration of simulated annealing with a cooling temperature 
ramp. This is the lowest level workhorse function used by both both flavors of simulated tempering that
have been developed.


#### use original "simulated tempering" scheduling algorithm
- *corrperm_ampdel_tempering.m* - function call from matlab to run M iterations of "simulated tempering" optimization
- *corrperm_ampdel_tempering_module.m* - compiled tempering module that can be run in a multiprocessing 
environment. This function reads its arguments from the command line provided by the MPE and passes them
to *corrperm_ampdel_tempering()*.
- *corrperm_tempering_schedule.m* - runs a chunk of multiple iterations 

#### use "lockstep" simulated annealing scheduling algorithm
This variant of the simulated tempering algorithm was used in [Schumacher2017] to better balance the fit between 
amplifications and deletions in the permutations with a small cohort of one cancer type. It allows multiple "reheats"
increase the better fitting of amps/dels more to favor the less better fitting alteration. 
- *corrperm_ampdel_lockstep.m* - function call from matlab to run M iterations of "lockstep" optimization
- *corrperm_ampdel_lockstep_module.m* - compiled "lockstep" module that can be run in multiprocessing 
environment. This function reads its arguments from the command line provided by the MPE and passes them
to *corrperm_ampdel_lockstep()*.
- *corrperm_lockstep_schedule.m* - runs a single iteration

#### adapters for different kinds of multiprocessing environments
- *corrperm_lsf_submission.m* - submit a module to (Broad Institute) Platform LSF
- *corrperm_lsf_wrapper.sh* - script to launch and pass parameters to compiled permutation module on one LSF CPU instance
- *corrperm_uger_submission.m* - submit a module to (Broad Institute) Grid Engine
- *corrperm_uger_wrapper.sh* - script to launch and pass parameters to compiled permutation module on one GridEngine CPU instance

#### data preparation
- *corrperm_prep.m* - prepare canonical reference data for permutations in a specified work directory from input copy number data in a D-struct and peaks. This main data preparation function was used in the analyisis for Zack2013.
- *corrperm_prep_gg.m* - variant of corrperm_prep() that prepares data for permutation considering the different background model for
deletions used by gene_gistic. Much more time-consuming than corrperm_prep.
- *create_D_sample_bins.m* - find the sample classes (cancer subtypes) in the input D-struct
- *score_D_scnas.m* - vestigial helper for corrperm_prep_gg()? Delete?
- *score_cooccurance.m* - create logical CN_event X sample arrays representing which sample has which events

### Analysis Functions
Analysis functions create a background distribution from the chunk file outputs of the controlled permutations
and assess the statistical significance of the observed data 
#### top-level functions
There are two top level analyses possible on permuted data: analyzing pairs of copy number events 
for anticorrelation or correlation, and analyzing the correlation (or anticorrelation) of a patient sample feature
with one other copy number event.
- *corrperm_analyze_pairs.m* - original algorithm for analyzing pairs of events which loads all permutations 
from all chunk files into memory at once. Fastest and most flexible for analysis development. Memory intensive: storage
scales NNP where *N* is the number of events and *P* is the number of permutations.
- *corrperm_analyze_pairs2.m* - more memory conservative analysis of event pairs, minimally scores permutations one chunk 
file at a time into a co-occurence count for each event. Storage scales independentky of number of permutations.
- *corrperm_analyze_pairs3.m* - slightly more flexible memory-conservative analysis, builds histogram of counts for each event from each chunk file, motivated primarily by the need to visualize.
- *corrperm_analyze_features.m* - original algorithm for analyzing features. This uses the same memory-intensive 
algorithm as the original *corrperm_analyze_pairs()* and should be revised similarly if memory is an issue in the analysis phase.
#### helper functions
- *max_fish_power.m* - calculate the maximum possible power (best p-value) for given marginals  
- *load_permutation_ref_inputs.m* - load the canonical reference data files into memory.
- *save_pair_pvalues.m* - write tab-delimited results file for correlated or anti-correlated pairs of CN events
- *save_feature_pvalues.m* - DELETE ME (replaced by the analyze_and_save() subfunction defined in corrperm_analyze_features.m)
- *save_feature_pvalues_2tailed.m* - DELETE ME (replaced by the analyze_and_save() subfunction defined in corrperm_analyze_features.m)
- *corrperm_get_final_stats.m* - this reads in the final amplification and deletion errors (difference from observed disruption) for each simulated annealing permutation from the permutation chunk files.

- *corrperm_display_stats.m* - displays the statistics gathered by the permutation.
