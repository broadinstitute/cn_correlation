function submit_perms(H,permuter,mpe,input_dir,perm_dir,Njobs,Niters,perm_opts)
%SUBMIT_PERMS submit permutation test chunks to multiprocessing environment
%
%    submit_perms(H,PERMUTER,MPE,INDIR,PERMDIR,NJOBS,NITERS,PERM_OPTS)
%
% PERMUTER name of the compiled matlab function that executes a chunk of permutations
% MPE names a file that defines the submit command for a multiprocessing environment 
% (e.g. lsf.sub or uger.sub)
% INPUT_DIR file path to directory containing inputs for permutation tasks
% PERM_DIR no longer used
% NJOBS is the number of chunks (schedulings of the permuter)
% NITERS is the number of permutations to perform per chunk
% PERM_OPTS is a struct of permutation options
%
% The inputs that need to be stored in REFDIR are:
%   - margs.mat - array of marginal disruption constraints for permutations
%   - new_samples.mat - the class definition for each sample
%
% The PERM_OPTS struct is stored as 'perm_ops.mat' in the PERM_DIR prior to
% submission.

if ~exist('perm_opts','var')
    perm_opts = struct;
    %!!! or should this be an error?
end

perm_opts.Nchunks = Njobs;
perm_opts.Niters = Niters;

perm_opts = impose_default_value(perm_opts,'output_subdir','permout'); %!!! move to function common to simulation and analysis
perm_opts = impose_default_value(perm_opts, 'verbose_level',40);

% make output directories as needed
if ~exist(input_dir,'dir')
    mkdir(input_dir);
end
% subdirectory for MPE outputs
output_dir = fullfile(input_dir,perm_opts.output_subdir);
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

% compile matlab executable and create matlabroot file
executable = corrperm_create_executable(input_dir,permuter);

% Find the wrapper script that sets up the matlab runtime environment 
wrapscript = which('matenvwrap.sh');
if isempty(wrapscript)
    error('cannot find matlab wrapper shell script ''matenvwrap.sh''');
end

% locate the multiprocessor environment specific file
platsub = which(mpe);
if isempty(platsub)
    error(['cannot find MPE platform command file ''',mpe,'''']);
end

mpe_cmd = which('mpe_cmd.sh');
if isempty(mpe_cmd)
    error('cannot find MPE variable setup file ''mpe_cmd.sh''');
end

if isfield(perm_opts,'randseed')
    rngseed = randseed(perm_opts.rngseed);
else
    rngseed = randseed();
    perm_opts.randseed = rngseed; 
end

% save permutation options in input_dir
save(fullfile(input_dir,'perm_opts.mat'),'perm_opts');
% save disruption profile
save(fullfile(input_dir,'H.mat'),'H');

%!save(fullfile(input_dir,'rngseed.mat'),'rngseed');

% create a launch script file in the reference directory
launch_script = fullfile(input_dir,'launch_perms.sh');
verbose('creating launch script ''%s''.',20,launch_script);

%% create launch script with an entry for each job
fid = fopen(launch_script,'w');
its = num2str(Niters);
for j = 1:Njobs
    chunk_id = sprintf('cycle%03d',j);  % chunk identifier extension string
    seed = num2str(randi(2^32-1));   % 32-bit chunk RNG seed string
    % Create a nested command to build a command to submit a job to a multi-processing environment using unix.
    % The first argument names a platform-generic shell script that sets some environment variables and then sources 
    % the next field, a platform-specific fragment that echos a "submit job" command for the MPE platform.
    % The third and remaining fields are the command and arguments to be run when the task is executed.
    % But the dance goes on: the first thing executed is a bash script that sets up the matlab environment,
    % then (finally!) launches the compiled matab executable along with its arguments.
    [s,cmdstr] = unix(strjoin({mpe_cmd,platsub,wrapscript,executable,input_dir,output_dir,chunk_id,its,seed}));
    if s ~= 0
        error('creating command from template failed');
    end
    fprintf(fid,cmdstr);
    verbose(cmdstr,30);
end
fclose(fid);

unix(['chmod u+x ',launch_script]);
%!!! optionally, defer launch so script can be executed manually
%!unix(launch_script);

end % function
