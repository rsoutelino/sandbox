
clear;clf;hold off;
close all

%*******************************************************
%                    PROGRAM MPSI_S
%
%          computes the Stommel's solution for specified
%            wind curl and bottom topography shapes
%
%            by Ilson da Silveira -- April 1997 
%
%            version with f-plane config added
%
%            by Ilson da Silveira -- August 2002 
%******************************************************

tol=1e-3;  % tolerance for convergence error

 
% parameters

L=1e6;
f0=7.3e-5;
beta=2e-11;
r=.5e-6/(beta*L); % bottom friction parameter
ff=f0/(beta*L);  

% set up grid

nmax=101;
dx=0.01; 

nsel1=[1:nmax-2];
nsel2=[3:nmax];
nsel=[2:nmax-1];


x=[0:nmax-1]*dx;
y=x';
[xg,yg]=meshgrid(x,y);
xg=flipud(xg);
yg=flipud(yg);

yL=max(y);


% wind curl matrix

%  curl=zeros(size(xg));
curl=-pi/yL*sin(pi*yg./yL);        % sinusoidal wind, single gyre
%  curl=pi/yL*sin(2*pi*yg./yL);     % sinusoidal wind, double gyre
%curl=-pi/yL*exp(-(xg-yL/2).*(xg-yL/2)); % gaussian
%curl(:,1:51)=-pi/yL*ones(101,51); % Heaviside step function
%curl=-pi/yL*ones(size(xg));

% nondimensional bottom topography--it is assumed H=5000m. 

   hb=zeros(size(xg)); % flat bottom
%  hb=-0.275*yg;      % linear sloped bottom in y;
%    hb=-3*0.275*xg;    % linear sloped bottom in x;
%  hb=-(xg-0.5);
%  hb=0.1*exp(-30*(xg-yL/2).*(xg-yL/2)); % gaussian mid-atlantic ridge
%  hb= 0.1*exp(-30*(xg-yL/2).*(xg-yL/2)).*exp(-30*(yg-yL/2).*(yg-yL/2));


% obtain gradients of ambient PV

fac=1.0                    % **set fac=0.0 for f-plane solution**
                           % otherwise, fac=1.0
 
y0=0.5;                    % reference latitude
B=fac*(yg-y0)+ff*hb;       % ambient PV function
%B=fac*(yg-y0)+hb;         % ambient PV function

[Bx,By]=gradient(B,dx,dx);  % zonal and meridional gradients of
By=-By;                     % ambient PV


gray2=jet(256);
%gray2=gray2(50:256-50,:);
colormap(gray2);

% plotting ambient PV

	    b1=0.1*floor(min(min(B))*10);
            b2=0.1*ceil(max(max(B))*10);


 lb=b1:.05:b2; 
 figure(1) 
 cs=contourf(x,y,flipud(B),lb,'w'); colorbar
% clabel(cs,lp)
 axis('square')
 xlabel('x')
 ylabel('y')
 title('The Ambient Potential Vorticity')




% plotting topography

% 2D view

 lt=-1:.05:1; 
 figure(2) 

 subplot(121)
 cs=contour(x,y,flipud(hb),lt,'k'); %colorbar
% clabel(cs,lp)
 axis('square')
 xlabel('x [in 10^6 m]')
 ylabel('y [in 10^6 m]')
% title('The Ambient Potential Vorticity')
 title('The Bottom Topography -- planar view')

% 3D view

 subplot(122)
 cs=surf(x,y,4e2*flipud(hb)); shading flat % colorbar
% clabel(cs,lp)
 axis('square')
 xlabel('x [in 10^6 m]')
 ylabel('y [in 10^6 m]')
 zlabel('z [in m]')
 title('The Bottom Topography -- 3D view')



% ITERATION PROCEDURE 

% boundary conditions

psi=zeros(size(xg));

% iteration begins

  for niter=1:2000

 niter
  
  psi_old=psi;
  av=psi(nsel1,nsel)+psi(nsel2,nsel)+psi(nsel,nsel1)+psi(nsel,nsel2);
  F=1/r*(curl(nsel,nsel)*dx*dx - psi(nsel,nsel2).*(By(nsel,nsel))*dx ...
  +psi(nsel2,nsel).*(Bx(nsel,nsel))*dx );
  psi(nsel,nsel)=(av-F)./(4+(dx/r)*(By(nsel,nsel)-Bx(nsel,nsel) ) );
  crit= max(max(abs(psi-psi_old)))
  if crit <= tol,break,end

  end

% compute relative vorticity

  zeta=(psi(nsel1,nsel)+psi(nsel2,nsel)+psi(nsel,nsel1)+psi(nsel,nsel2)- ...
  4*psi(nsel,nsel))./dx./dx;;


% plotting streamfunction

 orient portrait

	    p1=0.1*floor(min(min(psi))*10);
            p2=0.1*ceil(max(max(psi))*10);

 figure(3) 
 lp=p1:.1:p2; 
 cs=contourf(x,y,flipud(psi),lp,'w'); colorbar
% clabel(cs,lp)
 axis('square')
 xlabel('x [in 10^6 m]')
 ylabel('y [in 10^6 m]')
 title('The Stommel Solution')









