
function [tp,ep,A,C]=scaloa2(xc,yc,x,y,t,corrlen,err,zc)

%function [tp,ep]=scaloa2(xc,yc,x,y,t,corrlen,err,zc)
%	objectively maps the scalar variable t(x,y) to the points (xc,yc)
%	t can be m x n where columns hold values at (x,y) for particular
%	realization (e.g. day).
%       err is the normalized r.m.s.error (ICS, 05/06/93)
% 	If t is empty only the error values are produced.
%	If only tp is asked for no error values are produced.
%       version with two options for correlation fc types
%       remove bias a la Carter & Robinson '87
%       INPUT LINE VECTORS ONLY!

n=length(x);
x=reshape(x,1,n);
y=reshape(y,1,n);
t=t';

%array of squared distances between observation points

d2=((x(ones(n,1),:)'-x(ones(n,1),:)).^2+...
	(y(ones(n,1),:)'-y(ones(n,1),:)).^2);

nv=length(xc);
xc=reshape(xc,1,nv);
yc=reshape(yc,1,nv);

% array of squared distances between observation and interpolation points

dc2=((xc(ones(n,1),:)'-x(ones(nv,1),:)).^2+...
	(yc(ones(n,1),:)'-y(ones(nv,1),:)).^2);


if (nargin==7),

C=(1-err)*exp(-dc2/corrlen^2);
A=(1-err)*exp(-d2/corrlen^2);

else

C=(1-dc2./zc^2).*exp(-dc2/corrlen^2);
A=(1-d2./zc^2).*exp(-d2/corrlen^2);

end

A=A+err*eye(size(A));


if(~isempty(t))

         mD=sum(A\t)/sum(sum(inv(A)));
         t=t-mD;
	 tp=(C*(A\t))';
         tp=tp+mD*ones(size(tp));
end


if(nargout==4)
	ep=1-sum(C'.*(A\C'))/(1-err);
end





