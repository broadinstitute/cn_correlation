# cn_correlation

## Notes
- Original "simulated tempering" code for (Zack et al., *Nat Gen* 2013) was implemented on LSF multiprocessing. New 
"lockstep tempering" was implemented on GridEngine (Broad UGER) and used on (Schumacher et al., *PlosOne*, 2017) for 
smaller (600ish) data sets. 

## canonical reference files
- *Binary_amps.mat*
- *Binary_dels.mat*
- *new_samples.mat*
- *D.mat*
- *peak_regs.mat*

## Code File Descriptions
-README.txt - old documentation

### Permutation Functions

#### high level script examples
- *ng_pancan_corrperm_hl_recap.m* - script to recapitulate Zack "high-level" analysis using refactored functions
- *ng_pancan_corrperm_ll_recap.m* - script to recapitulate Zack "low-level" analysis using refactored functions

#### use "lockstep" simulated annealing scheduling algorithm
- *corrperm_ampdel_lockstep.m* - function call from matlab to run M iterations of "lockstep" optimization
- *corrperm_ampdel_lockstep_module.m* - compiled "lockstep" module that can be run in multiprocessing environment
- *corrperm_lockstep_schedule.m* - runs a chunk of multiple iterations of the "lockstep tempering"

#### use original "simulated tempering" scheduling algorithm
- *corrperm_ampdel_tempering.m* - function call from matlab to run M iterations of "simulated tempering" optimization
- *corrperm_ampdel_tempering_module.m* - compiled "simulated tempering" module that can be run in multiprocessing environment
- *corrperm_tempering_schedule.m* - runs a chunk of multiple iterations 

- *corrperm_cooldown_ratio.m* - run one iteration of simulated annealing with a temperature 
ramp. This is the lowest level workhorse function for both flavors of simulated tempering.
- *corrperm_get_final_stats.m*

#### adapters for different kinds of multiprocessing environments
- *corrperm_lsf_submission.m* - submit a module to (Broad Institute) Platform LSF
- *corrperm_lsf_wrapper.sh* - script to launch and pass parameters to compiled permutation module on one LSF CPU instance
- *corrperm_uger_submission.m* - submit a module to (Broad Institute) Grid Engine
- *corrperm_uger_wrapper.sh* - script to launch and pass parameters to compiled permutation module on one GridEngine CPU instance

#### data preparation helpers
- *corrperm_prep.m* - prepare canonical reference data for permutations in a specified work directory from input copy number data in a D-struct and peaks. This main data preparation function was used in the analyisis for Zack2013.
- *corrperm_prep_gg.m* - variant of corrperm_prep() that prepares data for permutation considering the different background model for
deletions used by gene_gistic. Much more time-consuming than corrperm_prep.
- *create_D_sample_bins.m* - find the sample classes (cancer subtypes) in the input D-struct
- *score_D_scnas.m* - vestigial helper for corrperm_prep_gg()? Delete?
- *score_cooccurance.m* - create logical CN_event X sample arrays representing which sample has which events

### Analysis Functions
#### top level functions
- *corrperm_analyze_features.m* -
- *corrperm_analyze_pairs.m* - original code for analyzing pairs of events which loads all permutations 
from all chunk files into memory at once. Fastest and most flexible alternative, but memory intensive.
- *corrperm_analyze_pairs2.m* - memory conservative analysis of event pairs, minimally scores permutations one chunk file at a time.
- *corrperm_analyze_pairs3.m* - slightly more flexible memory-conservative analysis, builds histogram of counts for each event.
#### helper functions
- *save_feature_pvalues.m* - write tab-delimited results file for correlation of CN events with some other patient feature
- *save_feature_pvalues_2tailed.m* -
- *save_pair_pvalues.m* - write tab-delimited results file for correlated or anti-correlated pairs of CN events
- *max_fish_power.m* - calculate the maximum possible power (best p-value) for given marginals  
- *load_permutation_ref_inputs.m* - load the canonical reference data files into memory.

- *corrperm_display_stats.m* - displays the statistics gathered by the permutation.
