
function [u,v] = psi2uv(x,y,psi)
% PSI2UV Velocity components from streamfunction
%   [U,V] = PSI2UV(X,Y,PSI) returns the velocity components U and V
%   from streamfunction PSI. X and Y are longitude and latitude matrices,
%   respectively, as PSI. All matrices have same dimension.
% 
%   --> !!! IMPORTANT !!!  <-----------------------------------------
%   The grid indexes IM and JM must have origin on the left-inferior 
%   corner, increasing to right and to up, respectively.
% 
%                              GRID
%            :            :            :            :
%            |            |            |            |
%   (JM = 2) o ---------- o ---------- o ---------- o --- ...
%            |            |            |            |
%            |            |            |            |
%   (JM = 1) o ---------- o ---------- o ---------- o --- ...
%            |            |            |            |
%            |            |            |            |
%   (JM = 0) o ---------- o ---------- o ---------- o --- ...
%        (IM = 0)     (IM = 1)     (IM = 2)     (IM = 3)
% 
%   -----------------------------------------------------------------
% 
%   USAGE:  [u,v] = psi2uv(x,y,psi)
% 
%   INPUT:
%     x     = decimal degrees (+ve E, -ve W) [-180..+180]
%     y     = decimal degrees (+ve N, -ve S) [- 90.. +90]
%     psi   = streamfunction [m^2 s^-1]
% 
%   OUTPUT:
%     u   = velocity zonal component [m s^-1]
%     v   = velocity meridional component [m s^-1]
% 
%   EXAMPLE:
%     [x,y] = meshgrid(-30:.5:-20,-35:.5:-25);
%     psi = (x-mean(mean(x))).^2 + (y-mean(mean(y))).^2;
%     [u,v] = psi2uv(x,y,psi);
%     contourf(x,y,psi,30);shading flat;hold on;
%     quiver(x,y,u,v,4,'w');

%   Author:
%   Rafael A. de Mattos (RAM), 03 May 2005
%   Update of 03 May 2005 (RAM)

%   ======================================================================


% Warning Off for regular grids
 
warning off MATLAB:divideByZero


[JM,IM] = size(psi);

u = zeros(JM,IM);
v = u;

% Velocity components u2 and v2 on Natural Coordinates

for i = 1:IM
  for j = 1:JM
  
    % Take for behind difference
    if j == 1
      dy2 = (y(j+1,i)-y(j,i))*60*1852;
      dx2 = (x(j+1,i)-x(j,i))*60*1852;
      dd = sqrt(dx2^2+dy2^2);
      V = (psi(j,i)-psi(j+1,i))/dd;
      if isinf(dy2/dx2) == 1
        ang = 0;
      elseif atan(dy2/dx2)*180/pi > 0 & isinf(dy2/dx2) == 0
        ang = atan(dy2/dx2)*180/pi - 90;
      else
        ang = atan(dy2/dx2)*180/pi + 90;
      end
      v2(j,i) = V*sin(ang*pi/180);
      u2(j,i) = V*cos(ang*pi/180);
    end
    
    
    % Take forward differences
    if j == JM
      dy2 = (y(j,i)-y(j-1,i))*60*1852;
      dx2 = (x(j,i)-x(j-1,i))*60*1852;
      dd = sqrt(dx2^2+dy2^2);
      V = (psi(j-1,i)-psi(j,i))/dd;
      if isinf(dy2/dx2) == 1
        ang = 0;
      elseif atan(dy2/dx2)*180/pi > 0 & isinf(dy2/dx2) == 0
        ang = atan(dy2/dx2)*180/pi - 90;
      else
        ang = atan(dy2/dx2)*180/pi + 90;
      end
      v2(j,i) = V*sin(ang*pi/180);
      u2(j,i) = V*cos(ang*pi/180);
    end
    
    % Take centered differences on interior points
    if j > 1 & j < JM
      dy2 = (y(j+1,i)-y(j-1,i))*60*1852;
      dx2 = (x(j+1,i)-x(j-1,i))*60*1852;
      dd = sqrt(dx2^2+dy2^2);
      V = (psi(j-1,i)-psi(j+1,i))/dd;
      if isinf(dy2/dx2) == 1
        ang = 0;
      elseif atan(dy2/dx2)*180/pi > 0 & isinf(dy2/dx2) == 0
        ang = atan(dy2/dx2)*180/pi - 90;
      else
        ang = atan(dy2/dx2)*180/pi + 90;
      end
      v2(j,i) = V*sin(ang*pi/180);
      u2(j,i) = V*cos(ang*pi/180);
    end
    
  end
end

% Velocity components u1 and v1 on Natural Coordinates

for j = 1:JM
  for i = 1:IM
  
    % Take for behind difference
    if i == 1
      dy1 = (y(j,i+1)-y(j,i))*60*1852;
      dx1 = (x(j,i+1)-x(j,i))*60*1852;
      dd = sqrt(dx1^2+dy1^2);
      V = (psi(j,i+1)-psi(j,i))/dd;
      ang = atan(dy1/dx1)*180/pi + 90;
      v1(j,i) = V*sin(ang*pi/180);
      u1(j,i) = V*cos(ang*pi/180);
    end
    
    
    % Take forward differences
    if i == IM
      dy1 = (y(j,i)-y(j,i-1))*60*1852;
      dx1 = (x(j,i)-x(j,i-1))*60*1852;
      dd = sqrt(dx1^2+dy1^2);
      V = (psi(j,i)-psi(j,i-1))/dd;
      ang = atan(dy1/dx1)*180/pi + 90;
      v1(j,i) = V*sin(ang*pi/180);
      u1(j,i) = V*cos(ang*pi/180);
    end
    
    % Take centered differences on interior points
    if i > 1 & i < IM
      dy1 = (y(j,i+1)-y(j,i-1))*60*1852;
      dx1 = (x(j,i+1)-x(j,i-1))*60*1852;
      dd = sqrt(dx1^2+dy1^2);
      V = (psi(j,i+1)-psi(j,i-1))/dd;
      ang = atan(dy1/dx1)*180/pi + 90;
      v1(j,i) = V*sin(ang*pi/180);
      u1(j,i) = V*cos(ang*pi/180);
    end
    
  end
end

% Components u and v on Cartesian Coordinates

v = v1+v2;
u = u1+u2;

return
