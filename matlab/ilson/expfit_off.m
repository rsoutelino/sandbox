function err = expfit(lambda)

%   err = expfit(lambda)
%   assumes an exponential function 
%
%     y =  c(1) + c(2)*exp(-lambda*t) 
%
%     c(1) is the fit offset
%
%   with the linear parameter vector c  and the  nonlinear parameter lambda.
%

global t y z err Plothandle

   A = exp(-lambda*t);
   A=[ones(size(A)) A];

c = A\y
z = A*c; 

set(Plothandle,'ydata',z)
drawnow
err = norm(z-y);
