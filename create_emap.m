function E = create_emap(D,regs,Binary_amps,Binary_dels);
%CREATE_EMAP create genomic disruption structure from processed GISTIC inputs
%
%    EMAP = create_emap(D,REGS,AMP_MATRIX,DEL_MATRIX)
%
% D is a GISTIC copy number structure
% REGS is a GISTIC peak region file
% AMP_MATRIX & DEL_MATRIX are logical calling matrix (Nevents x Nsamples)
%
% MARGS is a 3D array (Nchromosomes x Nsamples x amp/del type) of disruption scores
% PCX defines the permutation class structure of the data
% 
% The returned EMAP is an event mapping structure 
%
    % save event map for analysis as a single package
    E = struct;
    E.sdesc = D.sdesc;

    % event sub-structure
    E.event = struct;
    % take event names from peak names
    E.event.name = [strcat('amp_',{regs{1}.name}),strcat('del_',{regs{1}.name})]';
    % chromosome number from peaks
    E.event.chrn = [regs{1}.chrn,regs{2}.chrn]';
    % type: 1 = amp; 2 = del
    E.event.type = [repmat(1,length(regs{1}),1);repmat(2,length(regs{2}),1)];

    % event call matrix
    E.dat = [Binary_amps;Binary_dels];
 
end

