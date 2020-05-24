function seed = randseed(seed,varargin)
%RANDSEED initialize matlab random number generator for Monte Carlo
%
%     SEED = randseed(SEED,...)
%
% If called with no arguments or empty SEED, randomize using clock and return
% SEED used. Otherwise, initialize RNG using SEED. Additional arguments
% are passed to the rng() function (for Matlab version > 2010b).
%
% NOTE: matlab starts up with RNG in identical known state that must be 
%       deliberately randomized in order to use Monte Carlo techniques.

    % version-dependent behavior
    matlab8plus = str2num(regexprep(version,'\.[0-9]+\.[0-9]+ \(R.+\)$','')) >= 8.0;

    if ~exist('seed','var') || isempty(seed)
        % no seed: randomize random number generator for this chunk
        if matlab8plus
            rng('shuffle',varargin{:}); % Matlab R2012b and later
            seed = randi(2^32-1);
            rng(seed);
        else
            seed = rem(now*10e9,1e6); %! Matlab R2010b: no 'rng', seed rng from wall clock
            rand('seed',seed);
        end
    else
        % reproducible run: use passed-in seed
        if matlab8plus
            rng(seed,varargin{:});
        else
            rand('seed',seed);
        end
    end
end % function

