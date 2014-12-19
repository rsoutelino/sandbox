function [xi,yi,ui,vi]=mirror(xo,yo,uo,vo);

% function mirror computes the image of the
% velocity data points for OA streamline anlysis
%
% USAGE: [xi,yi,ui,vi] = mirror(xo,yo,uo,vo);
%
%     xo,yo are the corrdinates of the points
%
%     uo,yo are the zonal and meridional components of the velocity
%
%     xi,yi,ui,vi are the image points of the corrdinates and velocity components


% computing image position

m1=-(-45.5+53)/7.5;
b1= 7.5 -m1*(-53); 

m2=-1/m1;
b2= yo - m2*xo;

r= sqrt( m1*m1 + 1);

d= (-m1*xo + yo - b1)/r;

b3= -d*r +b1;

xi=-(b2-b3)/(m2-m1); 
yi= b3 + m1*xi;

%xi=xi-53;

% computing image velocity

th = atan(m2);

thv = atan(vo/uo);

vm=sqrt(uo*uo +vo*vo);

% 1st quadrant
if uo >=0 & vo >= 0,

   alf = th-thv;
   gam = pi + th + alf;

   ui=vm*cos(gam);
   
   vi=vm*sin(gam);

% 2nd quadrant

elseif uo < 0 & vo > 0, 
     thv = thv +pi;
    alf = abs(th - thv);
    gam = pi +th -alf;

      ui=vm*cos(gam);
   
      vi=vm*sin(gam); 

% 3rd quadrant

elseif uo <=0 & vo < 0,

   alf=th-thv;
   gam =2*pi +th+alf;

   ui=vm*cos(gam);
   
   vi=vm*sin(gam); 

% 4th quadrant
else
    th= th + pi;
    alf=abs(thv-th);
    gam=pi+th+alf;
 
   ui=vm*cos(gam);
   
   vi=vm*sin(gam); 
end






