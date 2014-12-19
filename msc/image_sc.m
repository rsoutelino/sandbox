function [xi,yi,psii] = image_sc(x,y,psi,xc,yc)
% IMAGE_SC Computes image points of scalar fields
%    [xi,yi,psii] = image_sc(x,y,psi,xc,yc) returns vectors
%    with image and real coordinates and values of scalar
%    fields (as streamfunction), to apply boudary conditions
%    to velocities using Objective Analisis interpolation
%    technique.
%
%  USAGE: [xi,yi,psii] = image_sc(x,y,psi,xc,yc)
%
%  INPUT:  
%           x = longitude vector, decimal degrees
%           y = latitude vector, decimal degrees
%         psi = scalar field vector
%          xc = x coordinates of the boundary (isobath)
%          yc = y coordinates of the boundary (isobath)
%           
%  OUTPUT: 
%          xi,yi = vectors containing real + image points
%          psii = vector containing real + image values
% 
%  IMPORTANT !!!
%         -  X and Y vectors need to be sorted as the first 
%  station of all transects is the one close to the coast. 
%         -  It only works for western boundaries! 
%  (the eastern, southern and northern boundaries will be 
%  implemented on the next version)
%
%  AUTHOR: 
%  Rafael G. Soutelino, 23 Dec 2006
%
%  =========================================================

dd = sw_dist(x,y,'km');
fd = find(dd > 200); % recognizing the transects

c = 0; lonest=[]; latest=[]; lonim=[]; latim=[]; psiest=[]; psiIM=[];
for k = 1:length(fd)+1
   
    if k == length(fd)+1;
        xrad = x(fd(end)+1:end);
        yrad = y(fd(end)+1:end);
        psirad = psi(fd(end)+1:end);
    else
        xrad = x((c+1):fd(k));
        yrad = y((c+1):fd(k));
        psirad = psi((c+1):fd(k));
        c = fd(k);
    end
    
    % cutting the real points that are outside the domain
    f = near(yc,yrad(1),1);
    f2 = find(xrad >= xc(f));
    xrad = xrad(f2); yrad = yrad(f2); psirad = psirad(f2);

    % colecting points of the boundary that are close to the transect
    xr = xc(f-2:f+2); yr = yc(f-2:f+2);
    p = polyfit(xr,yr,1); m1 = p(1); h1 = p(2); % y = m1.x + h1
    xr1 = min(x):0.1:max(x);
    yr1 = m1.*xr1 + h1;
    % imaging point by point
    xim=[]; yim=[]; psiim=[];
    for q = 1:length(xrad)
        m2 = -1/m1; % y = m2.x + h2
        h2 = xrad(q)/m1 + yrad(q);
        xi = (h2-h1)/(m1-m2); yi = m1*((h2-h1)/(m1-m2)) + h1; 
        u = xrad(q) - xi; v = yrad(q) - yi;
        x2 = xrad(q) - 2*u; y2 = yrad(q) - 2*v;

        % cutting image points that are inside the domain
        fn = near(yc,y2,1);
        if x2 < xc(fn) 
            xim = [xim x2]; yim = [yim y2]; psiim = [psiim -psirad(q)];
        end

    end 

   % mounting image and real vectors
   lonest = [lonest xrad]; latest = [latest yrad];
   lonim = [lonim xim]; latim = [latim yim];
   psiest = [psiest psirad]; psiIM = [psiIM psiim];

end

% mounting the final vectors:

xi = [lonest lonim];
yi = [latest latim];
psii = [psiest psiIM];

return