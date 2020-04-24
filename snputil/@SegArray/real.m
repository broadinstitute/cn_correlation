%    SegArray implementation of REAL.
%    B = REAL(A)
%    Emulates full array behavior of REAL; returns a SegArray.
%    (use HELP REAL for details)

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

function B = real(A)
% (code generated by unarray.pl)
B = A;
try
    B.vals = real(A.vals);
catch me
    throwAsCaller(me);
end