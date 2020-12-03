### clone the broadinstitute cn_correlation github repository
Recipe repositories are assumed to be under the directory ~/git. 
You should adjust the examples to match the location of your own
repository.

```
$ cd ~/git/broadinstitute
$ git clone https://github.com/broadinstitute/cn_correlation.git
$ cd cn_correlation/
```
### matlab setup
The user matlab startup script ~/matlab/startup.m runs when you start
up a matlab session. It is useful for setting up the paths matlab
searches for functions (defined in text files with the .m extension).
The recipe below assumes you are creating the file for the first
time: if you already have a startup file you need to edit it and add
the contents of the repo startup. Note that you should edit the file
anyway in order to reflect your own local git repo location.
```
$ cat startup.m
|if ~isdeployed
|    addpath ~/git/broadinstitute/cn_correlation
|    addpath ~/git/broadinstitute/cn_correlation/snputil
|end
$ mkdir ~/matlab
$ cp startup.m ~/matlab/
```
### copy the example directory files to a separate working directory (optional)
Although it is possible to work directly from the example directory, making a copy is
better because if you run a lot of permutations, the index storage can consume large 
amounts of disk space.
```
$ ls example/
$ mkdir /xchip/beroukhimlab/gistic/corrperm/example
$ cp example/* /xchip/beroukhimlab/gistic/corrperm/example
$ cd /xchip/beroukhimlab/gistic/corrperm/example
```

### customize the multiprocessing environment 'submit' template file
The submit template file creates the command for a multiprocessing environment to run
a task that will perform a "chunk" of permutations within the time frame the environment
allows for the task. A couple examples of submit template files are in the ```cn_correlation```
directory: lsf.submit is for LSF and has been tested at DFCI ErisOne; and one for GridEngine 
(UGER) tested at the Broad Institute.
```
$ cat ~/git/broadinstitute/cn_correlation/lsf.submit
| echo bsub -R "rusage[mem=8000]" -q medium -n 1 -o $OUTFILE -e $ERRFILE -r $*
$ cat ~/git/broadinstitute/cn_correlation/uger.submit
| echo qsub -N corrperm -b y -o $OUTFILE -e $ERRFILE -wd $WORK_DIR $*
```
These files may be edited, or copied and edited to suit your own multiprocessing environment. The following
environment strings are defined for the template:

- $* is the command to be executed
- $OUTPUT_DIR path to the output directory for permutations
- $WORK_DIR path to the input directory for permutations (soon this will be the same as $OUTPUT_DIR)
- $TAGOUT "output tag" that names the output associated with each specific chunk. For example if the output is named IDX_CELL.CYCLE001.mat the output tag is 'CYCLE001'
- $OUTFILE file to capture standard output from task
- $ERRFILE file to capture standard error from task

### do whatever you need to do on your cluster to run matlab 2014a
(Broad uses the dotkit 'use .matlab-2014a' command)

### then run matlab, and launch the example script
This reruns the 5000-permutation test on low-level disruption from
Travis's 2013 paper. It's worth looking at the first dozen and last
30 lines of ng_pancan_corrperm_ll_recap.m to understand what's going on.
The outputs will be placed in 3 directories
.../example/ll_work - input directory for the permutations
.../example/ll_work/ll_permout - output directory for the permutations
.../example/ll_work/ll_results - output directory for significance
                                 analysis of the permutations

```
$ matlab
>> ng_pancan_corrperm_ll_recap
```
If all goes well, this will run the data preparation and tuning phases for a few
minutes and then pause with a 'K>>' prompt from the 'keyboard'
matlab command. This allows you to enter matlab commands before
starting the permutations. If you type the command 'return' you will
run the edited corrperm_lsf_submission() function and attempt to
submit 100 permutation chunk tasks of 50 permutations each to LSF

K>> return

If by some miracle this works the first time, you will get another
keyboard 'K>>' prompt. Typing 'return' again will proceed to running 
the analysis corrperm_analyze_pairs2 once the completed permutations.
More likely you will keep editing corrperm_lsf_submission.m and
resubmitting until it works. To resubmit, type

>> corrperm_lsf_resub(ref_dir)

This matlab function picks up the parameters from the last LSF
submission stored in the file permute_options.mat in the reference 
directory and attempts a rerun. 

### analyze data

The last 10-ish lines of ng_pancan_corrperm_ll_recap do a "pairs2"
analysis of the data. I snipped them off and put them in
the file example/analyze_ll_pairs2.m to define a matlab command script you can
run by typing

>> analyze_ll_pairs2

                                                                                

