## cn_correlation inputs and outputs

### Permutation Module

#### Permutation Module Inputs

##### Main permutation inputs
This reflects the current state of things (2020) rather than the final plan.

*H* - disruption profile of input samples (TODO: higher, user-level description)
> The H struct contains a data matrix 'margs' of per-chromosome disruption per sample. The units are currently number of markers disrupted (problem: no total markers to normalize). The H struct also contains meta-data: 'sdesc' are sample_ids; 'chrname' are chromosome names; 'pcnames' and 'pcx' represent the lineage of each sample encoded as lineage names and corresponding indices of each lineage in 'sdesc' and 'margs'.

*permuter* - annealing/tempering algorithm to use, currently either *ampdel\_tempering\_module* or *ampdel\_lockstep\_module*

*mpe* - multi-processor environment customization file
>> currently either *lsf.submit* for the Platform LSF environment or *uger.submit* for the Univa GridEngine Research cluster

*ref_dir* - reference directory for input files
> Files written by permutation module:
>> *launch_perms.sh* - unix script to submit jobs to multiprocessor environment

>> *matlabroot* - text file containing path to matlab

>> compiled MCMC permuter module, currently either *ampdel\_tempering\_module* or *ampdel\_lockstep\_module*

*perm_dir* - permutation working/output directory
> Files written as annealing/tempering module inputs:
>> *perm_opts.mat* - MCMC module-specific options + 'randseed'

>> *H.mat* - cohort disruption profile

>> *rngseed* - seed for pseudo-random number generator

> Files written by compiled MCMC permuter module instances per-chunk:
>> *idx_cell.\<chunk>.mat* - array of sample reorderings calculated by the chunk

>> *stats.\<chunk>.mat* - MCMC progress statistics

>> *stat_finals.\<chunk>.mat* - final MCMC statistics for each chunk

> Files written by MPE capturing text output by each chunk:
>> *cycle\<chunk>.out.txt* - standard output

>> *cycle\<chunk>.err.txt* - standard error

*Njobs* - number of chunks to schedule with MPE

*Niters* - number of simulations per chunk

*params* - struct of input parameters, documented in next section

##### Permutation input parameters

*randseed* - seed used for pseudo-random number generator



#### Permutation Module Outputs

##### annealing/tempering module file inputs

*perm_opts.mat* - MCMC module-specific options + 'randseed'

*H.mat* - cohort disruption profile that serves as an ensemble to permute. This file is also an input to the analyses.

##### annealing/tempering module file outputs (per chunk)

*idx_cell.\<chunk>.mat* - arrays of sample reorderings

*stats.\<chunk>.mat* - MCMC progress statistics

*stat_finals.\<chunk>.mat* - final MCMC statistics for each chunk

##### MPE outputs (per chunk)

*cycle\<chunk>.out.txt* - standard output

*cycle\<chunk>.err.txt* - standard error

##### misc
*rngseed* - seed for pseudo-random number generator (reproducibility)



### Analysis Module
#### Types of Analyses
*pair_stats* - p- and q-value estimators for event pair correlation/anti-correlation
> Default analysis. Compares co-occurrences in the constrained permutations with the observed co-occurences to estimate p-values.

*feature_stats* - p- and q-value estimators for sample features vs events
> Sample features are defined by columns in the sample info table that, like events, have values 0 or 1. Note that this is an "add-on" analysis to *pair_stats*. It is activated by specifying the *sample_features* parameter to a column or list of columns in the sample_info table.

*ks_stats* - (PLANNED) p- and q-value estimators for sample quantities vs events
> Sample quantities are sample info columns with numeric values. This analysis will apply the Kolmogorov-Smirnov integrated difference statistic to determine if the distribution of samples with each event is different from those without. This will be an "add-on" analysis to *pair_stats*. It is activated by specifying the *sample_values* parameter to a column or list of columns in the sample_info table.

*histogram* - summarize co-occurrence counts across specified pairs
> Output for each pair is a table of counts of co-occurrences associated with the number of permutations where that count was achieved. Non-graphical output, an indicator column indicates if permutation count is 'eq','gt' or 'lt' observed. This analysis will be triggered by setting the *histogram_pairs* parameter to a list of event pairs of interest. Either event or both events of a pair can be set to \* to include all events, but the pair analysis is time-consuming. A feature or value name can be substituted for an event if *feature_stats* or *ks_stats* analyses are being run.
 
*subgroup* - create p-value estimators for each subset in a sample partition
> Introduced to breakdown the contributions of each individually permuted lineage to the p-value, this is actually more general and can be used with any categorical sample info column.

#### Analysis Module Inputs
##### Main analysis inputs
*event_map* - map of which samples (rows) have which events (columns)
> The first column is presumed to be the *sample\_id* key and its values should match the key in the *sample\_info* file to the The remaining columns 

*event_info* - information about each event
> Required columns: 'event' matching row naes of event map; 'chr' naming the chromosome the event is on; 'type', 'a' for amplifications, 'd' for deletions. 

*perm_dir* - input directory containing permutation module output files

*results_dir* - directory for analysis output files

*sample_info* - per-sample information
> One column is keyed into the event map sample column. The analysis is performed on 

*params* - struct of input parameters, documented in next section

##### Analysis input arameters

*ext* - extension for output files

*sig_thresh* - q-value cutoff for output files (default 0.25)

*power_thresh* - p-value threshold for power filter (default 0.1)
> The maximum power of an event pair is calculated from its marginal event frequencies using the Fisher exact test on the most extreme realiation of a contingency table. Pairs with insufficient power are excluded from the results before applying the Benjamini-Hochberg FDR calcxulation top obtain q-values.

*pcount* - pseudo-count used for pair p-value estimation (default 1)
> Added to numerator and denominator of p-value estimator = number of permuted event counts exceding observed count / total permutations.

*split_eq* - split how observed = permuted events boundary cases are counted (default false)
> I still find this practice questionable, it needs to be explained to me again.

*sample_features* - column(s) in the sample_info file defining features to test (default none)
> Setting this to a *sample_info* column or columns with logical values will trigger feature analyses of the co-occurence of the specified sample features with events.

*sample_values* - column(s) in the sample_info file defining values to test (default none)
> Setting this to a *sample_info* column or columns with numeric values will trigger Kolmogorov-Smirnov tests of each event to see if the distribution of values in samples with the event is different from the distribution of samples without. 

*sample_group* - column(s) in the sample_info file defining subgroups (default none)
> Setting this to a categorical *sample_info* column will trigger a subgroup "lineage" analysis for each value in the column.

*group_details* - enables group output details (default FALSE)
> If TRUE, then results for each subgroup as detailed as the overall analysis will be output, otherwise the output will be a table of p-values with dimensions event x group.

*do_pair_stats* - set to false to suppress event-pair analysis (default TRUE)

#### Analysis Module Outputs
