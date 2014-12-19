function chi = powerfit(coef0,data)
%  POWERFIT To fit a power function on data
%    CHI = POWERFIT(coef0,data) returns the sum difference squared 
%    between the curve fit and data. The power function has the form;
%
%                         F(X) = A.*X.^B + C
%
%    USAGE:  chi = powerfit(coef0,data)
%
%    INPUT:
%      coef0   = Coeficients [A B C]
%      data    = [X Y]
%
%    OUTPUT:
%      chi   =  Sum difference squared SUM((F(X)-Y).^2)
%

% Startparameters
A = coef0(1);
B = coef0(2);
C = coef0(3);

% Data (x,y)
X = data(:,1);
Y = data(:,2);

% Fit
FX = A.*X.^B+C;

% Difference: Data minus Fit
chi = sum((FX-Y).^2);


