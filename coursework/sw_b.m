
function beta = sw_b(lat)

% SW_F       BETA factor
%===========================================================================
% 
% USAGE:  beta = sw_b(lat)
%
% DESCRIPTION:
%    Calculates the Beta factor defined by
%       f = 2*Omega*cos(lat0)*y/a  where Omega = 7.292e-5 radians/sec
%                                    and a = 6378000 m (Earth radius)
%
% INPUT:  
%   lat = Central Latitude in decimal degress north [-90..+90]
%
%=========================================================================

% CALLER:  general purpose
% CALLEE:  none

%-------------
% CHECK INPUTS
%-------------
if nargin ~= 1
   error('sw_b.m:  Requires one input argument')
end %if  

%-------------
% BEGIN
%-------------
% Eqn p27.  Unesco 1983.
a=6378000;
DEG2RAD = pi/180;
OMEGA   = 7.292e-5;     %s-1   A.E.Gill p.597
beta    = 2*OMEGA*cos(lat*DEG2RAD)/a;

return
%===========================================================================

