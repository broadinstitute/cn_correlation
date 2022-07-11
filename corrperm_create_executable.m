function exepath = corrperm_create_executable(ref_dir,permuter)
%CREATE_EXECUTABLE recompile matlab module and update matlabroot if necessary
%
% corrperm_create_executable(ref_dir,permuter)
%
% Compile a matlab executable using current version of matlab if it is out of date,
% and write the matlab root directory.
%   REF_DIR is the path for the executable and matlabroot file outputs
%   PERMUTER is the name of executable module to compile, the path is searched for 
% a source file with extension '.m' added
%
% The matlabroot file stores the root directory of the current version of matlab
% and is used to set up the matlab runtime for a multiprocessing environment.


if isempty(ref_dir)
    error('empty reference directory');
end

if ref_dir(1) ~= filesep
    ref_dir = fullfile(pwd,ref_dir);
end

% make output directory if needed
if ~exist(ref_dir,'dir')
    mkdir(ref_dir);
end

% figure out if we need to recompile executable
need_compile = false;
% find source file
srcpath = which(permuter);
if isempty(srcpath)
    error('unable to find source file ''%s.m'' on path',permuter);
end

exepath = fullfile(ref_dir,permuter);
% determine root path for currently running version of matlab
[rv,mlr_text] = unix('which matlab | sed -e ''s/\/bin\/matlab//''');
if rv ~= 0
    error('unable to determine current version of matlab')
end

% figure out status of save matlab root file
mlr_file = fullfile(ref_dir,'matlabroot');
if ~exist(mlr_file,'file')
    % if it doesn't exist, recompile
    verbose('recompile needed: no matlabroot',20)
    need_compile = true;
else
    % if it doesn't match current version of matlab, recompile
    [rv,old_mlr] = unix(['cat ',mlr_file]);
    if rv ~= 0
        error('cannot access matlab root file ''%s''',mlr_file);
    end
    if ~strcmp(old_mlr,mlr_text)
        verbose('recompile needed: matlabroot changed',20)
        need_compile = true;
    end
end
% figure out status of executable
if ~exist(exepath)
    % if it doesn't exist, recompile
    verbose('recompile needed: no executable',20)
    need_compile = true;
else
    % if the source module is newer that the executable, recompile
    %@!!! NOTE: doesn't check date of other dependent source files)
    exeinfo = dir(exepath);
    srcinfo = dir(srcpath);
    if datenum(exeinfo.date) < datenum(srcinfo.date)
        verbose('recompile needed: source newer than executable',20)
        need_compile = true;
    end
end

% if we need to recompile, do so
if need_compile
    % save current directory and switch to output diretory
    save_dir = pwd;
    cd(ref_dir);
    try
        % compile the code
        mcc('-m',[permuter,'.m']); % (using full srcpath generates warning)
        % update matlab root file
        fid = fopen(mlr_file,'w');
        if fid > 0
            fprintf(fid,mlr_text);
            fclose(fid);
        else
            error('unable to write to matlabroot file ''%s''',mlr_file);
        end

        %!!!nope        unix(['echo ',mlr_text,'>',mlr_file]);
        % remove compiler-generated files we don't need
        unix('rm run_corrperm_ampdel_tempering_module.sh');
        unix('rm readme.txt');
        unix('rm mccExcludedFiles.log');
    catch me
        % switch back to original directory
        cd(save_dir);
        rethrow(me);
    end
    % no error: switch back to original directory
    cd(save_dir);
end
