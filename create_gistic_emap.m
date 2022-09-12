function E = create_gistic_emap(D,regs,pcx,Binary_amps,Binary_dels);
%CREATE_GISTIC_EMAP create genomic disruption structure from processed GISTIC inputs
%
%    EMAP = create_gistic_emap(D,REGS,AMP_MATRIX,DEL_MATRIX)
%
% inputs are outputs of corrperm_prep:
%   D is a GISTIC copy number structure
%   REGS is a GISTIC peak region file
%   AMP_MATRIX & DEL_MATRIX are logical matrices of event calls (Nevents x Nsamples)
%   PCX defines the permutation class structure of the data
% 
% The returned EMAP is an event mapping structure.
%
    % save event map for analysis as a single package
    E = struct;
    E.sample = struct('id',D.sdesc);

    % minimal event info
    E.event = struct('name',[strcat('amp_',{regs{1}.name}),strcat('del_',{regs{2}.name})]',...
                     'chrn',num2cell([regs{1}.chrn,regs{2}.chrn]'),...
                     'type',num2cell([repmat(1,length(regs{1}),1);repmat(2,length(regs{2}),1)] ));
    %{
    % take event names from peak names
    E.event.name = [strcat('amp_',{regs{1}.name}),strcat('del_',{regs{2}.name})]';
    % chromosome number and residual q-value from peaks
    E.event.chrn = [regs{1}.chrn,regs{2}.chrn]';
    E.event.resid_qv = [regs{1}.resid_qv,regs{2}.resid_qv]';
    % type: 1 = amp; 2 = del
    E.event.type = [repmat(1,length(regs{1}),1);repmat(2,length(regs{2}),1)];
    %}
    
    % lineage fields for analysis == permutation class fields (for now)
    %! E.pcx = pcx;
    %! E.pcname = extract_pcname(D,pcx);
    % event call matrix
    E.dat = [Binary_amps;Binary_dels];
 
end

