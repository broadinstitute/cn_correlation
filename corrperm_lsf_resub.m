function corrperm_lsf_resub(ref_dir)
%LSF_RESUBMIT resubmit permutation tests to LSF

% default reference directory is working directory
if ~exist('ref_dir','var')
    ref_dir = [pwd,'/'];
end

varfile = fullfile(ref_dir,'permute_options.mat');

if ~exist(ref_dir,'dir')
    error(['reference directory''',ref_dir,''' doesn''t exist!']);
end

if ~exist(varfile,'file')
    error(['''',varfile,''' does not exist!']);
end

% load resubmission parameters
load(varfile);
% opts,Njobs,Niters,permuter,perm_dir loaded


% resubmit
corrperm_lsf_submission(permuter,ref_dir,perm_dir,Njobs,Niters,opts);
