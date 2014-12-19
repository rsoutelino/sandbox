function coef=sinfitc(t,x,T,opt)

%**********************************************************************************
% function coef=sinfit(t,x,T,opt) fits a sinusoidal curve with 
% KNOWN PERIOD to a time series.
% 
% x(t)=dependent variable, t=time, T=selected period
%
% for opt=1, simple sinuosidal wave fit
%       newx(t)=coef(1)*sin((2*pi/T)*t)+coef(2)*cos((2*pi/T)*t) 
% two coefficientes are ouput
%
% for opt=2, simple sinuosidal wave fit + bias
%      newx(t)=coef(1)+coef(2)*sin((2*pi/T)*t)+coef(3)*cos((2*pi/T)*t)
% three coefficientes are ouput; the first is the bias
%
% for opt=3,  simple sinuosidal wave fit + bias +linear trend
%     newx(t)=coef(1)+coef(2)*t+coef(3)*sin((2*pi/T)*t)+coef(4)*cos((2*pi/T)*t) OR
% three coefficientes are ouput; the first is the bias; the second is the bias
%*************************************************************************************

xx=x(:);tt=t(:);
wt=2*pi/T;
%

if opt==1,
   A=[sin(wt*tt) cos(wt*tt)];
elseif opt==2
   A=[ones(size(tt)) sin(wt*tt) cos(wt*tt)];
else
A=[ones(size(tt)) ones(size(tt)).*tt sin(wt*tt) cos(wt*tt)];
end

coef=(A'*A)\(A'*xx);
