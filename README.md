# cn_correlation

## History
This algorithm was created by Travis Zack under the direction of Rameen Beroukhim. It was applied to a 5000 patient sampe pan-cancer cohort and published in [Zack213] (Zack et al., *Nat Gen*, 2013). Travis's working logs were adapted post-publication into a re-useable functional format by Steve Schumacher. The "lockstep tempering" variant was invented to improve the effectiveness of the algorithm with a smaller (700ish), single-disease cohort for [Schumacher2017] (Schumacher et al., *PlosOne*, 2017). The original implementation was adapted to use the GridEngine MPE (Broad UGER) and use less memory (allowing more permutations) for a 11K patient sample analysis published in [Gong2016] (Gong et al., *Neoplasia*, 2017). The algorithm has adapted for studies of correlation of arm-level changes in low grade gliomas (TCGA et al., *NEJM*, 2015) and aneuploidy (Taylor et al., *Cancer Cell*, 2018).

This code was originally maintained in the Broad Institute "CancerGenomeAnalysis" Subversion 1.6 repository under the path
https://svnrepos/CancerGenomeAnalysis/trunk/matlab/snp/correlation. It generally depends on some of the code under https://svnrepos/CancerGenomeAnalysis/trunk/matlab/

## Overview
-This method creates a background model for accessing the significance of observed co-occurring (or mutually excluding) copy number events by creating a large number of permutations from the observed data. 

## Workflow
- prepare data and measure amplification/deletion disruption for every patient sample's chromosome
- run a few tuning cycles
- run permutations in multi-processing environment
- analyze the permutations

### Core Algorithm

## Data Representation
Disruption is specified as the propensity to amplify or delete for each chromosome on each sample.

|disruption|sample1|sample2|...|sampleN|
|---|---|---|---|---|
|**chr1**| a<sub>11</sub> / d<sub>11</sub> | a<sub>12</sub> / d<sub>12</sub> |...| a<sub>1N</sub> / d<sub>1N</sub> |
|**chr2**| a<sub>21</sub> / d<sub>21</sub> | a<sub>22</sub> / d<sub>22</sub> |...| a<sub>2N</sub> / d<sub>2N</sub> |
|**chr3**| a<sub>31</sub> / d<sub>31</sub> | a<sub>32</sub> / d<sub>32</sub> |...| a<sub>3N</sub> / d<sub>3N</sub> |
|**chr4**| a<sub>41</sub> / d<sub>41</sub> | a<sub>42</sub> / d<sub>42</sub> |...| a<sub>4N</sub> / d<sub>4N</sub> |
|**chr5**| a<sub>51</sub> / d<sub>51</sub> | a<sub>52</sub> / d<sub>52</sub> |...| a<sub>5N</sub> / d<sub>5N</sub> |
|**chr6**| a<sub>61</sub> / d<sub>61</sub> | a<sub>62</sub> / d<sub>62</sub> |...| a<sub>6N</sub> / d<sub>6N</sub> |
|**...**|...|...|...|...|
|**chrX**| a<sub>X1</sub> / d<sub>X1</sub> | a<sub>X2</sub> / d<sub>X2</sub> |...| a<sub>XN</sub> / d<sub>XN</sub> |
|**genome**| A<sub>1</sub> / D<sub>1</sub> | A<sub>2</sub> / D<sub>2</sub> |...|A<sub>N</sub> / D<sub>N</sub> |

## permutation operations
- initial randomization
- objective function
- swap procedure

## implementation notes
- temperature units
- chunking
- subiterations

## canonical reference files output by corrperm_prep

### Analysis inputs

- *Binary_amps.mat* - sample x event logical matrix of amplification events
- *Binary_dels.mat* - sample x event logical matrix of deletion events
- *D.mat* - loaded, but not currently used
- *peak_regs.mat* - peak (event) definitions
- *new_samples.mat* - class (lineage) definition

### MPE permutation inputs
These are the identical inputs for each job.
- *margs.mat* - 3D array *margs_sort* of disruption values Nchr X Nsamples X 2 (amp/del) 
- *permute_options.mat* - contains a struct named *opts*
- *new_samples.mat* - class definition for the samples, cell array *new_samples* of index vectors

### MPE module per-chunk output files
Each MPE job outputs its chunk as four native matlab files.
- *rand_margs.{chunk}.mat* - randomized disruption matrices, currently not used by analysis
- *idx_cell.{chunk}.mat* - indices mapping original samples to permuted for each chromosome
- *stat_finals.{chunk}.mat* - final error for each run
- *stats.{chunk}.mat* - statistics for each run

## Code File Descriptions
-README.txt - old documentation, more thinking out how to functionally decompose Travis's scripts, less useful for final code.

### Permutation Functions

#### high level script examples
- *ng\_pancan\_corrperm\_hl\_recap.m* - script to recapitulate [Zack2013] "high-level" analysis using refactored functions
- *ng\_pancan\_corrperm\_ll\_recap.m* - script to recapitulate [Zack2013] "low-level" analysis using refactored functions

#### simulated annealing algorithm
- *corrperm\_cooldown\_ratio.m* - run one iteration of simulated annealing with a cooling temperature ramp. This is the lowest level workhorse function used by both both flavors of simulated tempering that have been developed.


#### use original "simulated tempering" scheduling algorithm
- *corrperm\_ampdel\_tempering.m* - function call from matlab to run M iterations of "simulated tempering" optimization
- *corrperm\_ampdel\_tempering\_module.m* - compiled tempering module that can be run in a multiprocessing 
environment. This function reads its arguments from the command line provided by the MPE and passes them
to *corrperm\_ampdel\_tempering()*.
- *corrperm\_tempering\_schedule.m* - runs a chunk of multiple iterations 

#### use "lockstep" simulated annealing scheduling algorithm
This variant of the simulated tempering algorithm was used in [Schumacher2017] to better balance the fit between amplifications and deletions in the permutations with a small cohort of one cancer type. It allows multiple "reheats"
increase the better fitting of amps/dels more to favor the less better fitting alteration. 
- *corrperm\_ampdel\_lockstep.m* - function call from matlab to run M iterations of "lockstep" optimization
- *corrperm\_ampdel\_lockstep\_module.m* - compiled "lockstep" module that can be run in multiprocessing 
environment. This function reads its arguments from the command line provided by the MPE and passes them
to *corrperm\_ampdel\_lockstep()*.
- *corrperm\_lockstep\_schedule.m* - runs a single iteration

#### adapters for different kinds of multiprocessing environments
- *corrperm\_lsf\_submission.m* - submit a module to (Broad Institute) Platform LSF
- *corrperm\_lsf\_wrapper.sh* - script to launch and pass parameters to compiled permutation module on one LSF CPU instance
- *corrperm\_uger\_submission.m* - submit a module to (Broad Institute) Grid Engine
- *corrperm\_uger\_wrapper.sh* - script to launch and pass parameters to compiled permutation module on one GridEngine CPU instance

#### data preparation
- *corrperm\_prep.m* - prepare canonical reference data for permutations in a specified work directory from input copy number data in a D-struct and peaks. This main data preparation function was used in the analyisis for Zack2013.
- *corrperm\_prep\_gg.m* - variant of corrperm_prep() that prepares data for permutation considering the different background model for
deletions used by gene_gistic. Much more time-consuming than corrperm_prep.
- *create\_D\_sample\_bins.m* - find the sample classes (cancer subtypes) in the input D-struct
- *score\_D\_scnas.m* - vestigial helper for corrperm_prep_gg()? Delete?
- *score\_cooccurance.m* - create logical CN_event X sample arrays representing which sample has which events

### Analysis Functions
Analysis functions create a background distribution from the chunk file outputs of the controlled permutations
and assess the statistical significance of the observed data 

#### top-level functions
There are two top level analyses possible on permuted data: analyzing pairs of copy number events 
for anticorrelation or correlation, and analyzing the correlation (or anticorrelation) of a patient sample feature
with one other copy number event.

- *corrperm\_analyze\_pairs.m* - original algorithm for analyzing pairs of events which loads all permutations 
from all chunk files into memory at once. Fastest and most flexible for analysis development. Memory intensive: storage
scales NNP where *N* is the number of events and *P* is the number of permutations.
- *corrperm\_analyze\_pairs2.m* - more memory conservative analysis of event pairs, minimally scores permutations one chunk 
file at a time into a co-occurence count for each event. Storage scales independentky of number of permutations.
- *corrperm\_analyze\_pairs3.m* - slightly more flexible memory-conservative analysis, builds histogram of counts for each event from each chunk file, motivated primarily by the need to visualize.
- *corrperm_analyze_features.m* - original algorithm for analyzing features. This uses the same memory-intensive 
algorithm as the original *corrperm_analyze_pairs()* and should be revised similarly if memory is an issue in the analysis phase.

#### helper functions
- *max\_fish\_power.m* - calculate the maximum possible power (best p-value) for given marginals  
- *load\_permutation\_ref\_inputs.m* - load the canonical reference data files into memory.
- *save\_pair\_pvalues.m* - write tab-delimited results file for correlated or anti-correlated pairs of CN events
- *save\_feature\_pvalues.m* - DELETE ME (replaced by the analyze\_and\_save() subfunction defined in corrperm\_analyze\_features.m)
- *save\_feature\_pvalues\_2tailed.m* - DELETE ME (replaced by the analyze\_and\_save() subfunction defined in corrperm\_analyze\_features.m)
- *corrperm\_get\_final\_stats.m* - this reads in the final amplification and deletion errors (difference from observed disruption) for each simulated annealing permutation from the permutation chunk files.

- *corrperm\_display\_stats.m* - displays the statistics gathered by the permutation.
