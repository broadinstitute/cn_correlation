function events = score_D_scnas(D,regs,t_amp,t_del,recon_params)
%SCORE_D_SCNAS - gene-gistic style event scoring
%
%   EVENTS = score_D_scnas(D,REGS,T_AMP,T_DEL,RECON_PARAMS)
%
% Returns a 1x2 cell array of amplification and deletion
% event matrices with a row for each peak specified in REGS.
% and a column for each sample in D, using the thresholds 
% T_AMP and T_DEL. Events are scored using the most extreme 
% value found in the boundaries specified by the 'peak_st' and
% 'peak_en' fields from the amp/del structs stored in the cell 
% array REGS. The D struct contains relative copy number (GISTIC
% intermediate result).
%
% If the RECON_PARAMS parameter is suplied, the scoring is performed on
% a genome reconstructed from D.Qs. The RECON_PARMS fields specifiy the
% focal reconstruction:
%   RECON_PARAMS.t_amp and RECON_PARAMS.t_del are noise thresholds
%   RECON_PARAMS.length_cutoff is the length cutoff. This can be
% different for amps versus dels by making it a cell array. If set > 2.0
% then all events will be used.
%
% note: the REGS input is the one stored in wide_peak_regs.mat,
% not the one from peak_regs.mat (why are they different?)


[Nmarkers,Nsamples] = size(D.dat);
amp_events = false(length(regs{1}),Nsamples);
del_events = false(length(regs{2}),Nsamples);
events = {amp_events,del_events};

% do not allow thresholds to be exactly zero (as they be for ABSOLUTE data)
t_amp = max(eps,t_amp);
t_del = max(eps,t_del);

thresh = [t_amp,-t_del];

% optional reconstruction
if exist('recon_params','var') && ~isempty(recon_params)
    % default noise thresholds are event calling thresholds
    if ~isfield(recon_params,'t_amp')
        recon_params.t_amp = t_amp;
    end
    if ~isfield(recon_params,'t_del')
        recon_params.t_del = t_del;
    end
    recon_params.use_segarray = true;
    recon_params.rows = Nmarkers;
    recon_params.cols = Nsamples;
    recon_params.broad_or_focal = 'focal';

    %% allow different length cutoffs for amps and dels
    if iscell(recon_params.broad_len_cutoff)
        %% reconstruct amps and dels with different length thresholds
        % amplification processing
        recon_params.broad_len_cutoff = recon_params.length_cutoff{1};
        focals = reconstruct_genomes(D.Qs,recon_params);
        events{1} = call_events(focals.amp,regs{1},t_amp);
        % deletion processing
        recon_params.broad_len_cutoff = recon_params.length_cutoff{2};
        focals = reconstruct_genomes(D.Qs,recon_params);
        events{2} = call_events(-focals.del,regs{2},-t_del);
    else
        % reconstruct with same amp and del length cutoff
        recon_params.broad_len_cutoff = recon_params.length_cutoff;
        focals = reconstruct_genomes(D.Qs,recon_params);
        focdat = focals.amp - focals.del;
        for k=1:2 % over amp/del
            events{k} = call_events(focdat,regs{k},thresh(k));
        end
    end
        
else

    %% call events without reconstruction
    for k=1:2   % loop over amp/del
        % loop over peaks in each SCNA type
        events{k} = call_events(D.dat,regs{k},thresh(k));
    end
end

%% low-level event calling function
%
% data_mat is marker X sample data matrix, peaks is a struct-array of peaks
% thresh is threshold (positive for amplifications, negative for deletions)
%
function event_mat = call_events(data_mat,peaks,thresh)
event_mat = false(length(peaks),size(data_mat,2));
for p = 1:length(peaks);
    markers = peaks(p).peak_st:peaks(p).peak_en;
    if thresh > 0
        event_mat(p,:) = max(data_mat(markers,:),[],1)>=thresh;
    else
        event_mat(p,:) = min(data_mat(markers,:),[],1)<=thresh;
    end
end

