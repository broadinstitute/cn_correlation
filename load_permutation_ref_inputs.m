function [D,Binary,Binary_amps,Binary_dels,new_samples,regs] = load_permutation_ref_inputs(input_folder)

load(fullfile(input_folder,'Binary_amps.mat'));
load(fullfile(input_folder,'Binary_dels.mat'));
Binary = [Binary_amps;Binary_dels];
load(fullfile(input_folder,'new_samples.mat')); % NOTE: contains sample reordering for each lineage
load(fullfile(input_folder,'D.mat'));
%!load(fullfile(input_folder,'samp_names.mat')); % NOTE: unused, may be useful for diagnostics
load(fullfile(input_folder,'peak_regs.mat')) %! new output
