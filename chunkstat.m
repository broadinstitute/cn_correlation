% chunk-processor base class defining general interface
% subclasses "plug in" to the analysis module
classdef chunkstat
    properties
        Nevents = 0;
        Nsamples = 0;
        Nperms = 0;
    end

    methods
        % constructor called before chunk processing
        % S = chunkstat(emap)
        % emap is an Nevents x Nsamples logical array of observed events in each sample
        function this = chunkstat(emap,varargin)
            [this.Nevents,this.Nsamples] = size(emap);
            verbose('%d events x %d samples observed',20,this.Nevents,this.Nsamples);
            this.Nperms = 0;
        end
        % called for every permutation
        function stat = perm(stat,emap)
            stat.Nperms = stat.Nperms + 1;
        end

        % called after chunk processing complete
        function done(stat,Npairs)
        end
    end
end
