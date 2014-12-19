function err=vmode(rd);

%*******************************************************************
% function err=vmode(rd);
%
% this function computes the QG pressure normal modes
%             for a given N2(z) profile
% the input required is the internal deformation radius (rd) in km
%
% the ouput 'err' is the residual for 'vmode' be zero, and therefore
% 1/rd^2 be an eigenvalue
%
% To find the rd value which zeroes the vmode function, do:
%
%   rd_true= fzero('vmode',rd_guess), 
%
%  where rd_guess is your initial guess for rd (in km)
%
%*******************************************************************

global n2f2;
global po;
global f;

n=length(n2f2);
d=zeros(1,n);
f=ones(1,n);

dz=(po(2)-po(1))*0.001;

nn=(n2f2(1:n-1)+n2f2(2:n))/2*dz;
dd=0;ff=1;g=1/rd/rd*dz;
for i=1:n-1
  dd=dd-g*ff;
  ff=ff+nn(i)*dd;
  d(i+1)=dd;
  f(i+1)=ff;
  end;
f=f(1:n-1);
f=f/sqrt(mean(f.*f));
xi0=mean(f.*f.*f);
%plot(f,-po);
err=dd;

