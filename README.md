# cn_correlation

## Notes
- Original "simulated tempering" code for (Zack et al., *Nat Gen* 2013) was implemented on LSF multiprocessing. New 
"lockstep tempering" was implemented on GridEngine (Broad UGER) and used on (Schumacher et al., *PlosOne*, 2017) for 
smaller (600ish) data sets. 

## File Descriptions
-README.txt - old documentation

### Permutation Functions

#### high level examples
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

- *corrperm_cooldown_ratio.m* - runs one iteration of simulated annealing with a temperature 
ramp. This is the lowest level workhorse function
- *corrperm_display_stats.m* - 
- *corrperm_get_final_stats.m*
- *corrperm_lsf_submission.m*
- *corrperm_lsf_wrapper.sh*

- *corrperm_prep.m*
- *corrperm_prep_gg.m*
- *create_D_sample_bins.m*

- *corrperm_uger_submission.m*
- *corrperm_uger_wrapper.sh*

- *farm_perms_lsf2_for_lineage_hl.m* - vestigial (DELETE ME)

- *load_permutation_ref_inputs.m*
- *score_D_scnas.m* - 
- *score_cooccurance.m* - 

###Analysis Functions
- *corrperm_analyze_features.m* -
- *corrperm_analyze_pairs.m* -
- *corrperm_analyze_pairs2.m* -
- *corrperm_analyze_pairs3.m* -

save_feature_pvalues.m
save_feature_pvalues_2tailed.m
save_pair_pvalues.m

max_fish_power.m

