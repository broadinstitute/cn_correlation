function [D,Binary,Binary_amps,Binary_dels,new_samples,regs] = grab_permutation_inputs(input_folder)

load([input_folder,'Binary_amps.mat']);
load([input_folder,'Binary_dels.mat']);
Binary = [Binary_amps;Binary_dels];
load([input_folder,'new_samples.mat']); % NOTE: contains sample reordering for each lineage
load([input_folder,'D.mat']);
%!load([input_folder,'samp_names.mat']); % NOTE: unused, may be useful for diagnostics
load([input_folder,'peak_regs.mat']) %! new output