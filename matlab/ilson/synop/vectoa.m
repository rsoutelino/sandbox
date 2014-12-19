
function psi=vectoa(xc,yc,x,y,u,v,corrlen,err,b)
%function psi=vectoa(xc,yc,x,y,u,v,corrlen,err)
%	(xc,yc) are vectors (row or column) of interpolation points
%	(x,y)   are vectors (row or column) of observation points
%	(u,v) 	are matrices of east and north components with each day
%		entered column-wise or row-wise
%	corrlen,err are correlation length scales and error for a
%		gaussian streamfunction covariance function
%	psi 	is streamfunction at	(xc,yc) returned in same
%		format as u or v
n=length(x);
x=reshape(x,1,n);% in case row vector
y=reshape(y,1,n);
[mu,nu]=size(u);
if(mu==n)
	u=[u;v];	%data entered column-wise
elseif(nu==n)
	u=[u';v'];	%data entered row-wise
else
	disp('Data array has incorrect size, does not match (x,y).')
	return
end
% angles and distances
t=atan2(-y(ones(n,1),:)'+y(ones(n,1),:),-x(ones(n,1),:)'+x(ones(n,1),:));
d2=((x(ones(n,1),:)'-x(ones(n,1),:)).^2+...
	(y(ones(n,1),:)'-y(ones(n,1),:)).^2);
lambda=1/corrlen^2;
bmo=b*err/lambda;
R=exp(-lambda*d2);	%longitudinal
S=R.*(1-2*lambda*d2)+bmo;	%transverse
R=R+bmo;

A=zeros(2*n,2*n);
A(1:n,1:n)=(cos(t).^2).*(R-S)+S;
A(1:n,n+1:2*n)=cos(t).*sin(t).*(R-S);
A(n+1:2*n,1:n)=A(1:n,n+1:2*n);
A(n+1:2*n,n+1:2*n)=(sin(t).^2).*(R-S)+S;
A=A+err*eye(size(A));

% angles and distances
[nv1,nv2]=size(xc);
nv=nv1*nv2;
xc=reshape(xc,1,nv);%in case row vector
yc=reshape(yc,1,nv);
tc=atan2(-yc(ones(n,1),:)'+y(ones(nv,1),:),-xc(ones(n,1),:)'+x(ones(nv,1),:));
d2=((xc(ones(n,1),:)'-x(ones(nv,1),:)).^2+...
	(yc(ones(n,1),:)'-y(ones(nv,1),:)).^2);
R=exp(-lambda*d2)+bmo;

P=zeros(nv,2*n);	%streamfunction-velocity covariance
P(:,1:n)=sin(tc).*sqrt(d2).*R;
P(:,n+1:2*n)=-cos(tc).*sqrt(d2).*R;

psi=(P*(A\u))'; %uses this line if u is full
if(nu==n)		%adjust if column-wise
	psi=psi';
end
