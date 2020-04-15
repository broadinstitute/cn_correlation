function [ Binary,Binary_amps,Binary_dels ] = score_cooccurance(Qs,amps,dels,params,do_broad,hl)
%SCORE_COOCCURANCE - call amp/del events at peak regions in copy-nmber data
% [EVENTS,AMP_EVENTS,DEL_EVENTS] = SCORE_COOCCURANCE(Qs,amps,dels,params,do_broad,hl)
%
% Qs is a matrix of events (see Qs)
% AMPS is a vector of amplification SNP indices
% DELS is a vector of deletion SNP indices
% PARAMS is a struc of parameter values for reconstruction
% DO_BROAD is a flag, if 0 purely focal events are used for scoring, if 1,
%          all events above the threshold are used for scoring
% HL is the SCNA threshold (above noise level) used to call events

%   Detailed explanation goes here

%fields = {'amp','del','aod','doa'};    
%Qs_cur = rmfield(Qs,fields);

%% amplification events
% reconstruct and apply threshold
focal_genome = reconstruct_genomes(Qs,params);
F.dat = focal_genome.amp-focal_genome.del;
Binary_amps = F.dat(amps,:)>(params.t_amp+hl);
%% deletion events
if do_broad
    % reconstruct using all events above noise threshold for deletions
    params.broad_len_cutoff = 2.1;
    broad_genome = reconstruct_genomes(Qs,params);
    B.dat = broad_genome.amp-broad_genome.del;
    Binary_dels = B.dat(dels,:)<-(params.t_del+hl);
else
    % threshold focal events for deletions like amplifications
    Binary_dels = F.dat(dels,:)<-params.t_del;
end
Binary_amps = full(Binary_amps);
Binary_dels = full(Binary_dels);
Binary = [Binary_amps;Binary_dels];
end

