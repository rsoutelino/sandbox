function w = window(N,wt)
%
%  w = window(N,wt)
%
%  generate a window function
%
%  N = length of desired window
%  wt = window type desired
%       'rect' = rectangular        'tria' = triangular (Bartlett)
%       'hann' = Hanning            'hamm'  = Hamming
%       'blac' = Blackman
%
%  w = row vector containing samples of the desired window
nn = N-1;
pn = 2*pi*(0:nn)/nn;
if wt(1,1:4) == 'rect',
                        w = ones(1,N);
elseif wt(1,1:4) == 'tria',
                        m = nn/2;
                        w = (0:m)/m;
                        w = [w w(ceil(m):-1:1)];
elseif wt(1,1:4) == 'hann',
                        w = 0.5*(1 - cos(pn));
elseif wt(1,1:4) == 'hamm',
                        w = .54 - .46*cos(pn);
elseif wt(1,1:4) == 'blac',
                        w = .42 -.5*cos(pn) + .08*cos(2*pn);
else
                        disp('Incorrect Window type requested')
end
