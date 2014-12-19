function [c,mu,sigma] = gaussfit_example(data)
%Function is an example how to use the gaussfit.m Function
%  data:  Data [x; y]

% set start parameters: c mu sigma
param0 = [1 1 1];

% call fminsearch
par = fminsearch(@gaussfit,...
               param0,...
               optimset('display','on','tolx',1e-6,'tolfun',1e-6),...
               data,1);

% read best parameters	       
c = par(1,1);           
mu = par(1,2);
sigma = par(1,3);

disp(['    c: ' num2str(c)]);
disp(['   mu: ' num2str(mu)]);
disp(['sigma: ' num2str(sigma)]);
