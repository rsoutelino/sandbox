function xf = gauss_filter2D(x,lh,lv,n)
% GAUSS_FILTER2D  2-D non-dimensional gaussian filter
%    XF = GAUSS_FILTER2D(X,LH,LV,N) filters the field contained in
%    the matrix X applying a weighed mean in the center point. The
%    weights are distributed as a bidimensional gaussian shape with
%    horizontal width LH and vertical width LV. This function assumes
%    the matrix domain is non-dimensional, e.g., it ranges between
%    0 and 1. To N = 1 the weight distribution is ploted to each 
%    point of the domain.
%

%   Author:
%   Rafael A. de Mattos (RAM), 26 February 2007
%   Update of 26 February 2007 (RAM)

%   ======================================================================

[l,c] = size(x);

[xg,yg] = meshgrid(0:1/(c-1):1,1:-1/(l-1):0);

for I = 1:c
  for J = 1:l
    C = exp(-((xg-xg(J,I)).^2./(lh).^2+(yg-yg(J,I)).^2./(lv).^2));
    if n == 1
      figure(1);contourf(xg,yg,C,20);shading flat;
    end
    xf(J,I) = sum(sum(x.*C))./sum(sum(C));
  end
end