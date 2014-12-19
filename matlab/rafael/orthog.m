
function [ort,xo,yo] = orthog(x,y)
% ORTHOG Orthogonality of the grid
%   [ORT,XO,YO] = ORTHOG(X,Y) returns the orthogonality of the grid 
%   where X and Y are longitude and latitude matrices, respectively. 
%   The orthogonality is represented by the deviation degree from 
%   orthogonality on the intersection point between four grid cell, 
%   exclusive of vertices along the edge of the grid.
% 
%               : y
%          /    :     \
%         o     :      o
%        / \    :     / \
%       /   \   :    /   \
%      /     \  ;   /     \            - orthogonality:
%             \ ;  /
%              \; / a1                   |a1| + |a2| = 90 [degree]
%               ;/\...............
%               ;\/             x      - function return:
%              /; \ a2
%             / ;  \                     ORT = 90 - (|a1| + |a2|)
%      \     /  ;   \     /
%       \   /        \   /
%        \ /          \ /
%         o            o
%          \          /
% 
% 
%   USAGE:  [ort,xo,yo] = orthog(x,y)
% 
%   INPUT:
%     x     = decimal degrees (+ve E, -ve W) [-180..+180]
%     y     = decimal degrees (+ve N, -ve S) [- 90.. +90]
% 
%   OUTPUT:
%     ort   = decimal degrees (deviation) [-90..90]
%     xo    = decimal degrees (+ve E, -ve W) [-180..+180]
%     yo    = decimal degrees (+ve N, -ve S) [- 90.. +90]

%   Author:
%   Rafael A. de Mattos (RAM), 03 May 2005
%   Update of 03 May 2005 (RAM)

%   ======================================================================

% Warning Off for regular grids
 
warning off MATLAB:divideByZero

[JM,IM] = size(x);

ort = NaN*ones(JM,IM);

for i = 2:IM-1
  for j = 2:JM-1
    
    dy2 = y(j,i+1)-y(j,i-1);
    dx2 = x(j,i+1)-x(j,i-1);
    ang2 = atan(dy2/dx2)*180/pi;
    
    dy1 = y(j+1,i)-y(j-1,i);
    dx1 = x(j+1,i)-x(j-1,i);
    ang1 = atan(dy1/dx1)*180/pi;
    
    ort(j,i) = 90-(abs(ang2)+abs(ang1));
        
  end
end

x(isnan(ort)) = [];xo = reshape(x,JM-2,IM-2);
y(isnan(ort)) = [];yo = reshape(y,JM-2,IM-2);
ort(isnan(ort)) = [];ort = reshape(ort,JM-2,IM-2);

return
