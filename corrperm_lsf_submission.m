function corrperm_lsf_submission(permuter,ref_dir,perm_dir,Njobs,Niters,opts)
%!!! what will be function inputs
%!ref_dir = '/xchip/gistic/schum/mcmc_correlation_test_130930/work_dir/';
%!perm_dir = '/xchip/gistic/schum/mcmc_correlation_test_130930/output_dir/';
%!permuter = 'corrperm_ampdel_tempering_module';
%!Njobs = 10; % number of jobs to run
%!Niters = 20; % number of iterations

% wrapper script !!!TODO existence test
wrapscript = which('corrperm_lsf_wrapper.sh');

srcpath = which(permuter);
if isempty(srcpath)
    throw(MException('matlab:submit_perm_jobs:no_src',...
                     'cannot find source file ''%s''.',permuter));
else
    % check for executable in pwd
    executable = fullfile(pwd,permuter);
    if exist(executable,'file')
        % check to see if recompile is needed ('make' functionality)
        exeinfo = dir(executable);
        srcinfo = dir(srcpath);
        if datenum(exeinfo.date) < datenum(srcinfo.date)
            mcc('-m',srcpath);
            %! TODO!!! use deptree to get all dependencies for date
            %checking
        end
    else
        % executable does not exist - compile it
        mcc('-m','-v',srcpath);
    end
end

cmd_template = ['bsub ',...
          '-R "rusage[mem=4]" ',... 
          '-mig 5 ',...
          '-R "select[cpuf>100]" ',... 
          '-Q "EXCLUDE(127)" ',...
          '-q hour ',...
          '-W 240 ',...
          '-P cancerfolk ',...
          '-o ',perm_dir,'%s.out.txt -e ',perm_dir,'%s.err.txt ',...
          '-r ',wrapscript,' ',executable,' ',ref_dir,' ',perm_dir ' %s %d']; 
      
      
% need input directory
if ~exist(ref_dir,'dir')
    error('input directory doen''t exist!')
end
save(fullfile(ref_dir,'permute_options.mat'),'opts');
% make output directory if need be
if ~exist(perm_dir,'dir')
    mkdir(perm_dir);
end
% submit to the farm
for j=1:Njobs
    runid = sprintf('cycle%03d',j);
    unix_str = sprintf(cmd_template,runid,runid,runid,Niters);
    disp(unix_str);
    [r1,r2]=unix(unix_str); 
    disp([strtrim(r2) ' Exit code: ' num2str(r1)])
end