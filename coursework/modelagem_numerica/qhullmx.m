function [t,v] = qhullmx(X, flags)
% QHULLMX gateway function for Qhull.
%   T = QHULLMX(X, 'd ', 'QJ') this is equivalent to T = DELAUNAYN(X).
%
%   [T,V] = QHULLMX(X, 'v ', 'QJ', 'o') is equivalent to [T,V] = VORONOIN(X).
%
%   QHULLMX is based on Qhull.  See QHULL for Copyright information.
%
%   See also QHULL, DELAUNAYN, VORONOIN.

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.7 $ $Date: 2002/06/05 20:52:32 $
%# mex  

error('mex-file not found')
