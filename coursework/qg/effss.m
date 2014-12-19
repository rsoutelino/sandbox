function nb=effss(x);

% function nv=effss(x);
% estimates effective sample size
% or # of degrees of freedom of x

x=x(:);
n=length(x);

if n<2,
  disp('error: length(x) < 2')
  return
end

  rho=sum( x(2:n).*x(1:(n-1)) )/sum( x(2:n).^2 );
nb=n*((1-rho)/(1+rho));
%nv=n*((1-rho^2)/(1+rho^2));
