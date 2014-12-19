function chi = gaussfit(param0,data,boolPlot)
%Function fits a Gauss-Function to Data (x,y)
%  y=c*exp(-(x-mu)^2/(2*sigma^2))
%Its parameters are: c,mu,sigma.
%
%function chi = gaussfit(param,data,boolPlot)
%   param0: Start parameters c,mu,sigma
%   data:   Data [x; y]
%   boolPlot:  1=realtime plotting,
%              0=no realtime plotting
%

% Startparameters
c = param0(1);
mu = param0(2);
sigma = param0(3);

% Data (x,y)
x = data(1,:);
y = data(2,:);

% Fit: Gauss-Function
fx = c.*exp(-((x-mu).^2)./(2*sigma.^2));

% Difference: Data minus Fit
chi = sum((fx-y).^2); 

% Realtime plotting
if boolPlot==1
    % plot Data
    plot(x,y,'b');
    hold on
    % plot Fit
    plot(x,fx,'r');
    hold off
    drawnow
end
