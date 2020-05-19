function corrperm_lsf_submission(permuter,ref_dir,perm_dir,Njobs,Niters,opts)
%LSF_SUBMISSION submit permutation test to Platform Load Sharing Facility
%
% lsf_submission(permuter,refdir,permdir,Njobs,Niters,options)
%
% PERMUTER name of the matlab function that executes a chunk of permutations
% REFDIR file path to directory containing inputs
% PERMDIR file path to directory where output chunk files will be stored
% NJOBS is the number of chunks (schedulings of the permuter)
% NITERS is the number of permutations to perform per chunk
% OPTS is a struct containing PERMUTER-dependent tuning parameters
%
% The inputs that need to be stored in REFDIR are:
%   - margs.mat - array of marginal disruption constraints for permutations
%   - new_samples.mat - the class definition for each sample
%
% The OPTS struct is stored as 'permute_options.mat' in the REFDIR prior to
% submission, along with a 'launch_perms' shell script.
 
% need input directory
if ~exist(ref_dir,'dir')
    error('input (reference) directory doen''t exist!')
end

% make output directory if need be
if ~exist(perm_dir,'dir')
    mkdir(perm_dir);
end

% save annealing parameters in opts, plus other stuff for ease of resubmission
save(fullfile(ref_dir,'permute_options.mat'),'opts','permuter','perm_dir','Njobs','Niters');

% wrapper script !!!TODO existence test
wrapscript = which('corrperm_lsf_wrapper.sh');
if isempty(wrapscript)
    error('cannot find wrapper shell script for LSF');
end

% compile matlab executable and create matlabroot file
executable = corrperm_create_executable(ref_dir,permuter);

cmd_template = ['bsub ',...
          '-R "rusage[mem=4]" ',... 
          '-mig 5 ',...
          '-R "select[cpuf>100]" ',... 
          '-Q "EXCLUDE(127)" ',...
          '-q SHORT ',...
          '-W 240 ',...
          '-P cancerfolk ',...
          '-o ',perm_dir,'/%s.out.txt -e ',perm_dir,'/%s.err.txt ',...
          '-r ',wrapscript,' ',executable,' ',ref_dir,' ',perm_dir ' %s %d\n']; 
      
launch_script = [ref_dir,'launch_perms.sh'];
verbose('creating launch script ''%s''.',20,launch_script);
fid = fopen(launch_script,'w');
for j=1:Njobs
    runid = sprintf('cycle%03d',j);
    fprintf(fid,cmd_template,runid,runid,runid,Niters);
end
fclose(fid);
unix(['chmod u+x ',launch_script]);

% bsub chunk tasks to LSF
for j=1:Njobs
    runid = sprintf('cycle%03d',j);
    unix_str = sprintf(cmd_template,runid,runid,runid,Niters);
    disp(unix_str);
    [r1,r2]=unix(unix_str); 
    disp([strtrim(r2) ' Exit code: ' num2str(r1)])
end

