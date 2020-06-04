function submit_perms(H,permuter,mpe,ref_dir,perm_dir,Njobs,Niters,perm_opts)
%SUBMIT_PERMS submit permutation test chunks to multiprocessing environment
%
%    submit_perms(H,PERMUTER,MPE,REFDIR,PERMDIR,NJOBS,NITERS,PERM_OPTS)
%
% PERMUTER name of the compiled matlab function that executes a chunk of permutations
% MPE names a file that defines the submit command for a multiprocessing environment 
% (e.g. lsf.sub or uger.sub)
% REF_DIR file path to directory containing inputs
% PERM_DIR file path to directory where output chunk files will be stored
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

% need input directory
if ~exist(ref_dir,'dir')
    error('input directory doen''t exist!')
end

% make output directory if need be
if ~exist(perm_dir,'dir')
    mkdir(perm_dir);
end

% save permutation options in perm_dir
save(fullfile(perm_dir,'perm_opts.mat'),'perm_opts');
save(fullfile(perm_dir,'H.mat'),'H');

% compile matlab executable and create matlabroot file
executable = corrperm_create_executable(ref_dir,permuter);

% Find the wrapper script that sets up the matlab runtime environment 
wrapscript = which('matenvwrap.sh');
if isempty(wrapscript)
    error('cannot find matlab wrapper shell script ''matenvwrap.sh''');
end

platsub_file = 'lsf.submit'; %!!! from settings
platsub = which(mpe);
if isempty(platsub)
    error(['cannot find MPE platform command file ''',mpe_plat_file,'''']);
end

mpe_cmd = which('mpe_cmd.sh');
if isempty(mpe_cmd)
    error('cannot find MPE variable setup file ''mpe_cmd.sh''');
end

if isfield(perm_opts,'randseed')
    rngseed = randseed(perm_opts.rngseed);
else
    rngseed = randseed();
end
save(fullfile(perm_dir,'rngseed.mat'),'rngseed');

% create a launch script file in the reference directory
launch_script = fullfile(perm_dir,'launch_perms.sh');
verbose('creating launch script ''%s''.',20,launch_script);

%% create launch script with an entry for each job
fid = fopen(launch_script,'w');
its = num2str(Niters);
for j = 1:Njobs
    chunk_id = sprintf('cycle%03d',j);  % chunk identifier extension string
    seed = num2str(randi(2^32-1));   % 32-bit chunk RNG seed string
    % Create a nested command to build a command to submit a job to a multi-processing environment using unix.
    % The first command names a platform-generic shell script that sets some environment variables and then sources 
    % the next field, a platform-specific fragment that echos a "submit job" command for the MPE platform.
    % The third and remaining fields are the command and arguments to be run when the task is executed.
    % But the dance goes on: the first thing executed is a bash script that sets up the matlab environment,
    % then (finally!) launches the compiled matab executable along with its arguments.
    [s,cmdstr] = unix(strjoin({mpe_cmd,platsub,wrapscript,executable,ref_dir,perm_dir,chunk_id,its,seed}));
    if s ~= 0
        error('creating command from template failed');
    end
    fprintf(fid,cmdstr);
end
fclose(fid);

unix(['chmod u+x ',launch_script]);
%!!! optionally, defer launch so script can be executed manually
unix(launch_script);

end % function
