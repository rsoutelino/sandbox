function [m,d,h,min] = decjul2greg(y,t)

%==================================================================================
% decjul2greg.m - Converts the decimal days of CODAS database into gregorian days
%
%  USAGE: 
%           [m,d,h,min] = decjul2greg(y,t)
%
%  INPUT: 
%           t   =    time in decimal days  (1 x 1  or  1 x m)
%           y   =    year of the database  (1 x 1)  (e.g. 2007) 
%
%  OUTPUT: 
%           m   =    month     (1-12)
%           d   =    day       (1-31)
%           h   =    hour      (1-23)  
%         min   =    minute    (1-59)
%
%  by:                Rafael Soutelino
%  last update:       28/08/2007
%==================================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin > 2
  error('decimal2greg.m: No more than 1 argument allowed')
end

% checking if the year has 365 or 366
 
if rem(y,4)==0 
   bix = 1;
   mon = [1:12 ; 31 29 31 30 31 30 31 31 30 31 30 31];
   fmon = [cumsum(mon(2,:))];
else
   bix = 0;
   mon = [1:12 ; 31 28 31 30 31 30 31 31 30 31 30 31];
   fmon = [cumsum(mon(2,:))];
end

lt = length(t);
M=[]; D=[]; H=[]; MIN=[];

for k = 1:length(t)

tk = t(k);

% finding the month
f = near(fmon,tk,1);

if fmon(f) < tk
   m = f + 1;
   fd = fmon(f);
else
   m = f;
   fd = fmon(f-1);
end  

M = [M m];

% finding the day
rd = floor(tk);
d = rd - fd+1;

D = [D d];

% finding the hour
h = tk-rd;
hd = h*24;
h = floor(hd);

H = [H h];

% finding the minutes
min = hd-h;
min = floor(min*60);

MIN = [MIN min];

clear min hd h tk d rd fd m f

end

clear min hd h tk d rd fd m f

m=M; d=D; h=H; min=MIN;

return







