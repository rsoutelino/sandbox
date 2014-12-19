function hpol = mypolar2(theta,rho,line_style,fillcolour,rinci,rmaxi)
%MYPOLAR2	Polar coordinate plot.
%	POLAR(THETA, RHO) makes a plot using polar coordinates of
%	the angle THETA, in radians, versus the radius RHO.
%	POLAR(THETA,RHO,S,fc) uses the linestyle specified in string S.
%       fc is a string which specifies fill colour - if there is
%       only 2 input arguments no fill is done.
%       rinci is the radial increment, rmaxi is the maximum radial
%       distance from the origin.  rinci and rmaxi are optional, and
%       rinci can be put in without rmaxi. PM 23 Sept 1996.
%       At present the fill overwrites the circles and 
%       ylabel.  
%	See PLOT for a description of legal linestyles.
%
%	See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.

%	Copyright (c) 1984-94 by The MathWorks, Inc.
%% PM changed the orientation;
%% to give 0 as the top (to match with north).
%% and increasing theta in a clockwise direction.
theta = -theta+pi/2;

if nargin < 1
	error('Requires 2 or 4 input arguments.')
elseif nargin == 2 
	if isstr(rho)
		line_style = rho;
		rho = theta;
		[mr,nr] = size(rho);
		if mr == 1
			theta = 1:nr;
		else
			th = (1:mr)';
			theta = th(:,ones(1,nr));
		end
	else
		line_style = 'auto';
	end
elseif nargin == 1
	line_style = 'auto';
	rho = theta;
	[mr,nr] = size(rho);
	if mr == 1
		theta = 1:nr;
	else
		th = (1:mr)';
		theta = th(:,ones(1,nr));
	end
end
if isstr(theta) | isstr(rho)
	error('Input arguments must be numeric.');
end
if any(size(theta) ~= size(rho))
	error('THETA and RHO must be the same size.');
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
	'DefaultTextFontName',   get(cax, 'FontName'), ...
	'DefaultTextFontSize',   get(cax, 'FontSize'), ...
	'DefaultTextFontWeight', get(cax, 'FontWeight') )

% only do grids if hold is off
if ~hold_state

% make a radial grid
	hold on;
	hhh=plot([0 max(theta(:))],[0 max(abs(rho(:)))]);
	v = [get(cax,'xlim') get(cax,'ylim')];
	ticks = length(get(cax,'ytick'));
	delete(hhh);
% check radial limits and ticks
	rmin = 0; rmax = v(4); rticks = ticks-1;
         if nargin ==6
          %%PM added this to have the optional input variable
          % to change the rmax. (6th variable in the function call).
          rmax = rmaxi;
         end

	if rticks > 5	% see if we can reduce the number
		if rem(rticks,2) == 0
			rticks = rticks/2;
		elseif rem(rticks,3) == 0
			rticks = rticks/3;
		end
	end

% define a circle
	th = 0:pi/50:2*pi;
	xunit = cos(th);
	yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = [1:(length(th)-1)/4:length(th)];
    xunits(inds(2:2:4)) = zeros(2,1);
    yunits(inds(1:2:5)) = zeros(3,1);

	rinc = (rmax-rmin)/rticks;
        if nargin==5| nargin==6
        % PM added this to put in optional input of the
        % rinc variable (rinci is the 5th variable in the
        % function call).
         rinc = rinci;
        end

	for i=(rmin+rinc):rinc:rmax
		plot(xunit*i,yunit*i,'-','color',tc,'linewidth',1);
		text(0,i+rinc/20,['  ' num2str(i)],'verticalalignment','bottom' );
	end

% plot spokes
	th = (1:6)*2*pi/12;
	cst = cos(th); snt = sin(th);
	cs = [-cst; cst];
	sn = [-snt; snt];
	plot(rmax*cs,rmax*sn,'-','color',tc,'linewidth',1);

% annotate spokes in degrees-PM changed to annotating in N,S,W,E.
 	rt = 1.1*rmax;
	for i = 1:max(size(th))
%%PM		text(rt*cst(i),rt*snt(i),int2str(i*30),'horizontalalignment','center' );
          if i==3 
           text(rt*cst(i),rt*snt(i),'N','horizontalalignment','center' );
           text(-rt*cst(i),-rt*snt(i),'S','horizontalalignment','center' );
          end
          if i==6 
           text(rt*cst(i),rt*snt(i),'W','horizontalalignment','center' );
           text(-rt*cst(i),-rt*snt(i),'E','horizontalalignment','center' );
          end

		if i == max(size(th))
			loc = int2str(0);
		else
			loc = int2str(180+i*30);
		end
%%PM		text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center' );
	end
% set viewto 2-D
	view(0,90);
% set axis limits
	axis(rmax*[-1 1 -1.1 1.1]);
end 

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
	'DefaultTextFontName',   fName , ...
	'DefaultTextFontSize',   fSize, ...
	'DefaultTextFontWeight', fWeight );

% transform data to Cartesian coordinates.
xx = rho.*cos(theta);
yy = rho.*sin(theta);
% plot data on top of grid

  %PM added the next section to fill the polar plot with
  % colour.
 %if nargin==4|nargin==5|nargin==6
 % for ip=0:length(xx)-3
 %  part1 = [0 xx(ip+2) xx(ip+3) 0];
 %  part2 = [0 yy(ip+2) yy(ip+3) 0];
 %  fill(part1,part2,fillcolour)
 % end
%
% end
	rinc = (rmax-rmin)/rticks;
        if nargin==5| nargin==6
        % PM added this to put in optional input of the
        % rinc variable (rinci is the 5th variable in the
        % function call).
         rinc = rinci;
        end

	for i=(rmin+rinc):rinc:rmax
		plot(xunit*i,yunit*i,'-','color',tc,'linewidth',1);
		text(0,i+rinc/20,['  ' num2str(i)],'verticalalignment','bottom' );
	end
  % end of section added by PM 12 June 1996.

if strcmp(line_style,'auto')
	q = plot(xx,yy);
else
	q = plot(xx,yy,line_style);
end
if nargout > 0
	hpol = q;
end
if ~hold_state
	axis('equal');axis('off');
end

% reset hold state
if ~hold_state, set(cax,'NextPlot',next); end





