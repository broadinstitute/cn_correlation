%    SegArray implementation of REALLOG.
%    B = REALLOG(A)
%    Emulates full array behavior of REALLOG; returns a SegArray.
%    (use HELP REALLOG for details)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function B = reallog(A)
% (code generated by unarray.pl)
B = A;
try
    B.vals = reallog(A.vals);
catch me
    throwAsCaller(me);
end