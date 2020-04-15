function [a_err,d_err] = corrperm_get_final_stats(filespec)
% 17-Sep-2014 gather run statistics

files = dir(filespec);
nfiles = length(files);
a_err = [];
d_err = [];
n = 0;
for i = 1:nfiles
    x = load(['reference/perms/',files(i).name]);
    if ~isfield(x,'stat_finals') || ~iscell(x.stat_finals)
        warning('unrecognized stat_finals file format');
        break;
    end
    if ~exist('ad_err','var')
        pperf = length(x.stat_finals);
        ad_err = nan(pperf*nfiles,2);
        n = 0;
    end
    ad_err(n+1:n+pperf,:) = vertcat(x.stat_finals{:});
    n = n+pperf;
end

if exist('ad_err','var')
    a_err = ad_err(:,1);
    d_err = ad_err(:,2);
end