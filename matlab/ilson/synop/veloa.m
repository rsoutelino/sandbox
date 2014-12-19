function [vc,ec]=veloa(xc,yc,x,y,u,v,corrlen,err)
% 
%	(xc,yc) are vectors (row or column) of interpolation points
%	(x,y)   are vectors (row or column) of observation points
%	(u,v) 	are matrices of east and north components with each day
%		entered column-wise or row-wise
%	corrlen,err are correlation length scales and error for a
%		gaussian streamfunction covariance function
%       vc      has two columns containing the interpolated non-divergent
%               velcities
%       ec      is the normalized r.m.s. of the non-divergent velocity 
%               magnitude
%
n=length(x);
x=reshape(x,1,n);% in case row vector
y=reshape(y,1,n);
[mu,nu]=size(u);
if(mu==n)
	u=[u;v];	%data entered column-wise
elseif(nu==n)
	u=[u';v'];	%data entered row-wise
	v=v';
else
	disp('Data array has incorrect size, does not match (x,y).')
	return
end

% angles and distances
t=atan2(-y(ones(n,1),:)'+y(ones(n,1),:),-x(ones(n,1),:)'+x(ones(n,1),:));
d2=((x(ones(n,1),:)'-x(ones(n,1),:)).^2+...
	(y(ones(n,1),:)'-y(ones(n,1),:)).^2);
lambda=1/corrlen^2;
R=exp(-lambda*d2);	%longitudinal
S=R.*(1-2*lambda*d2);	%transverse

A=zeros(2*n,2*n);
A(1:n,1:n)=(cos(t).^2).*(R-S)+S;
A(1:n,n+1:2*n)=cos(t).*sin(t).*(R-S);
A(n+1:2*n,1:n)=A(1:n,n+1:2*n);
A(n+1:2*n,n+1:2*n)=(sin(t).^2).*(R-S)+S;
A=A+err*eye(A);

% angles and distances
[nv1,nv2]=size(xc);
nv=nv1*nv2;
xc=reshape(xc,1,nv);%in case row vector
yc=reshape(yc,1,nv);
tc=atan2(-yc(ones(n,1),:)'+y(ones(nv,1),:),-xc(ones(n,1),:)'+x(ones(nv,1),:));
d2=((xc(ones(n,1),:)'-x(ones(nv,1),:)).^2+...
	(yc(ones(n,1),:)'-y(ones(nv,1),:)).^2);
R=exp(-lambda*d2);
S=R.*(1-2*lambda*d2);	%transverse


C=zeros(2*nv,2*n);
C(1:nv,1:n)=(cos(tc).^2).*(R-S)+S;
C(1:nv,n+1:2*n)=cos(tc).*sin(tc).*(R-S);
C(nv+1:2*nv,1:n)=cos(tc).*sin(tc).*(R-S);
C(nv+1:2*nv,n+1:2*n)=(sin(tc).^2).*(R-S)+S;

td=atan2(-yc(ones(nv,1),:)'+yc(ones(nv,1),:),-xc(ones(nv,1),:)'+xc(ones(nv,1),:));
d2c=((xc(ones(nv,1),:)'-xc(ones(nv,1),:)).^2+...
	(yc(ones(nv,1),:)'-yc(ones(nv,1),:)).^2);

R=exp(-lambda*d2c);
S=R.*(1-2*lambda*d2c);	%transverse
B=zeros(2*nv);
B(1:nv,1:nv)=((cos(td).^2).*(R-S)+S);
B(1:nv,nv+1:2*nv)=cos(td).*sin(td).*(R-S);
B(nv+1:2*nv,1:nv)=B(1:nv,nv+1:2*nv); 
B(nv+1:2*nv,nv+1:2*nv)=(sin(td).^2).*(R-S)+S;
B=B + err*eye(B);


vc=(C*(A\u))';

ec=sqrt(  1 - ( C*(A\C'))./B );
ec=reshape(ec,nv*nv,1);


if(nu==n)		%adjust if column-wise
	vc=vc';
end






