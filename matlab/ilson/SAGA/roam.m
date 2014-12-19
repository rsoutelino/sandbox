function  aa=roam(a,incx,incy)
%  Roam(a,incx,0) moves the matrix 'a' 'incx' elements
%	to the right (left if <0), zero-padding the empty elements and so 
%	keeping the same matrix size.
%  Roam(a,0,incy) moves the matrix 'a' 'incy' elements
%	down (up if <0), zero-padding the empty elements and so 
%	keeping the same matrix size.
[ny,nx]=size(a);
%
if abs(incx)>nx|abs(incy)>ny,
 error('[incx,incy] must be less than size(a).')
end  
% Move in +x direction
aa=a;
%
if incx>0,
 aa=[zeros(ny,incx) a(1:ny,1:nx-incx)];
elseif incx<0,
 aa=[a(1:ny,1-incx:nx) zeros(ny,-incx)];
end
%
if incy>0,
 aa=[zeros(incy,nx); aa(1:ny-incy,1:nx)];
elseif incy<0,
 aa=[aa(1-incy:ny,1:nx); zeros(-incy,nx)];
end



