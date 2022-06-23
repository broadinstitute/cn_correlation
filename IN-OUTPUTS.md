## cn_correlation inputs and outputs

### Permutation Module

### Analysis Module
#### Types of Analysis
*pair_stats* - p- and q-value estimators for event pair correlation/anti-correlation
> Default analysis. Compares co-occurrences in the constrained permutations with the observed co-occurences to estimate p-values.

*feature_stats* - p- and q-value estimators for sample features vs events
> Sample features are defined by columns in the sample info table that, like events, have the value 0 or 1.

*ks_stats* - (PLANNED) p- and q-value estimators for sample quantities vs events
> Sample quantities are sample info columns with numeric values. This analysis will apply the Kolmogorov-Smirnov integrated difference statistic to determine if the distribution of samples with each event is different from those without.

*histogram* - summarize co-occurrence counts across specified pairs
> Output for each pair is a table of counts of co-occurrences associated with the number of permutations where that count was achieved. Non-graphical output, an indicator column indicates if permutation count is 'eq','gt' or 'lt' observed.
 
*lineage* - break down the 

#### Analysis Inputs
##### Files
*event_map* - map of which samples (rows) have which events (columns)
> The first column is presumed to be the sample_id key and its values should match the key in the sample_info file to the The remaining columns 

*event_info* - information about each event
> Required columns: 'event' matching row naes of event map; 'chr' naming the chromosome the event is on; 'type', 'a' for amplifications, 'd' for deletions. 

*perm_dir* - directory containing permutation module output files

*results_dir* - directory for analysis output files

*sample_info* - per-sample information
> One column is keyed into the event map sample column. The analysis is performed on 

##### Parameters

*ext* - extension for output files

*sig_thresh* - q-value cutoff for output files (default 0.25)

*power_thresh* - p-value threshold for power filter (default 0.1)
> The maximum power of an event pair is calculated from its marginal event frequencies using the Fisher exact test on the most extreme realiation of a contingency table. Pairs with insufficient power are excluded from the results before applying the Benjamini-Hochberg FDR calcxulation top obtain q-values.

*pcount* - pseudo-count used for pair p-value estimation (default 1)
> Added to numerator and denominator of p-value estimator = number of permuted event counts exceding observed count / total permutations.

*split_eq* - split how observed = permuted events boundary cases are counted
> I still find this practice questionable, it needs to be explained to me again.
