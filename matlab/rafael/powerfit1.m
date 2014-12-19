function chi = powerfit(coef0,data)
%  POWERFIT To fit a power function on data
%    CHI = POWERFIT(coef0,data) returns the sum difference squared 
%    between the curve fit and data. The power function has the form;
%
%                           F(X) = A.*X.^B
%
%    USAGE:  chi = powerfit(coef0,data)
%
%    INPUT:
%      coef0   = Coeficients [A B]
%      data    = [X Y]
%
%    OUTPUT:
%      chi   =  Sum difference squared SUM((F(X)-Y).^2)
%

% Startparameters
A = coef0(1);
B = coef0(2);

% Data (x,y)
X = data(:,1);
Y = data(:,2);

% Fit
FX = A.*X.^B;

% Difference: Data minus Fit
chi = sum((FX-Y).^2);


