function corrperm_uger_submission(permuter,ref_dir,perm_dir,Njobs,Niters,opts)
%UGER_SUBMISSION submit permutation test to Univa Grid Engine for Research
%
% uger_submission(permuter,refdir,permdir,Njobs,Niters,options)
%
% PERMUTER name of the matlab function that executes a chunk of permutations
% REF_DIR file path to directory containing inputs
% PERM_DIR file path to directory where output chunk files will be stored
% NJOBS is the number of chunks (schedulings of the permuter)
% NITERS is the number of permutations to perform per chunk
% OPTS is a struct containing PERMUTER-dependent tuning parameters
%
% The inputs that need to be stored in REFDIR are:
%   - margs.mat - array of marginal disruption constraints for permutations
%   - new_samples.mat - the class definition for each sample
%
% The OPTS struct is stored as 'permute_options.mat' in the REFDIR prior to
% submission.

% need input directory
if ~exist(ref_dir,'dir')
    error('input directory doen''t exist!')
end

% save annealing parameters in opts, plus other stuff for ease of resubmission
save(fullfile(ref_dir,'permute_options.mat'),'opts','permuter','perm_dir','Njobs','Niters');

% make output directory if need be
if ~exist(perm_dir,'dir')
    mkdir(perm_dir);
end

% compile matlab executable and create matlabroot file
executable = corrperm_create_executable(ref_dir,permuter);

% wrapper script
wrapscript = which('corrperm_uger_wrapper.sh');
if isempty(wrapscript)
    error('cannot find wrapper shell script for UGER');
end

% create a launch script file in the reference directory
launch_script = fullfile(ref_dir,'launch_perms.sh');
verbose('creating launch script ''%s''.',20,launch_script);
fid = fopen(launch_script,'w');
fprintf(fid,'qsub -t 1-%d -wd %s %s %s %s %s %d\n',...
        Njobs,perm_dir,wrapscript,executable,ref_dir,perm_dir,Niters);
fclose(fid);
unix(['chmod u+x ',launch_script]);


% submit to GridEngine as an array job (-t 1-N)
unix_str = ['qsub -t 1-',num2str(Njobs),' -wd ',perm_dir,...
            ' ',wrapscript,' ',executable,' ',ref_dir,' ',perm_dir,' ',num2str(Niters)];

fprintf('submitting qsub command:\n  ''%s''\n',unix_str);
[r1,r2]=unix(unix_str);
disp([strtrim(r2) ' Exit code: ' num2str(r1)]);
