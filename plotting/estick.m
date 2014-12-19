function h = estick(tt,x,y,tlen,units,labels)
% ESTICK        modified STICKPLOT from the oceans toolbox
%
% h = estick(tt,x,y,tlen,units,labels)
%
% STICKPLOT This function will produce a stickplot of vector data in 
%           chronological order. By calling STICKPLOT(DATE,VX,VY,TLEN), 
%           vectors with components (VX,VY) at times DATE are arranged in 
%           rows of length TLEN. Everything is scaled "nicely" and a 
%           scale-bar is drawn, with units given by an  (optional 5th 
%           parameter) UNITS string. 
%
%           If TLEN is a 2-element vector, then only that range of 
%           days is plotted.
%
%           If VX,VY are matrices, then column will be plotted above one
%           another (i.e. for time series at different depths). In this
%           case TLEN *must* be a 2-element vector. Finally, a 6th
%           parameter LABELS will be used to label each line.
% 
% Modifications:
%
%
% See also FEATHER QUIVER (although these are really useless because of
% their inability to get the right direction of arrows unless AXIS EQUAL)

% This is a really kludgey piece of code - my contribution to
% global confusion! 
%              -RP 
%
% Yes, but it works really good. A few sqrews needed tightening,
% though. All modifications are put inside ---- lines or marked "%%% e".
%                                                       Even
%
%Time-stamp:<Last updated on 03/06/02 at 13:31:56 by even@gfi.uib.no>
%File:</home/janeven/matlab/evenmat/estick.m>

[N,M]=size(x);

% ---- check inarguments (even) ----------------------
error(nargchk(3,6,nargin));
if nargin < 6 | isempty(labels)
  labels=[];
end
if nargin < 5 | isempty(units)
  units='units';
end
if nargin < 4 | isempty(tlen)
  tlen=mima(tt); 
end
% ----------------------------------------------------

if (min(N,M)==1),
  stackseries=0;
  if (nargin<6), M=0; end;
else
  stackseries=1;
  if (nargin<6),
    labels=[];
    if (M>1),
      for i=1:M,
        labels=strvcat(labels,int2str(i));   %%% e
        % strvcat is more robust than []-concatenation %%% e
      end;
    else
      labels='';
    end;
  end;
end;

% convert to column vectors
tt=tt(:);
x=x(:);
y=y(:);

%%if (nargin<4),             %%% e (replaced above)
%%  tlen=max(tt)-min(tt);
%%end;
%%if (nargin<5),
%%  units='units';
%%end;

if (max(size(tlen))==1),
   if (stackseries), error('TLEN must be a 2-element vector'); end;
   autoaxis=1;
   leftx=min(tt);
else
   autoaxis=0;
   leftx=tlen(1);
   tlen=tlen(2)-tlen(1);
end;



maxmag=max(max(sqrt(x.*x+y.*y)));
sc=tlen/maxmag/8;
wx=x*sc;            % scaled versions of the vectors
wy=y*sc;

% make this stuff stackable in rows of length tlen
if (autoaxis), xax=tlen+2*max(wx);
else           xax=tlen; end;

yax=xax;
if (autoaxis),
   yoff=yax/ceil( (tt( length(tt))-tt(1))/tlen);

   t=tt-floor( (tt-tt(1))/tlen)*tlen;
   yy=-yoff/2-floor( (tt-tt(1))/tlen)*yoff;
else
   if (stackseries),
      t=tt*ones(1,M);
      t=t(:);
      yoff=yax/M;
      yy=-yoff/2-ones(size(tt))*[0:M-1]*yoff;
      yy=yy(:);
   else
      yoff=yax;
      t=tt;
      yy=-yoff/2*ones(size(tt));
   end;
end;

xp=[ t  t+wx  NaN*ones(size(t)) ]';xp=xp(:);
yp=[ yy wy+yy NaN*ones(size(t)) ]';yp=yp(:);

% Now plot
%%clf reset;
h=plot(xp,yp,'-');
%h=line(xp,yp,'linestyle','-');

maxmag=10^round(log10(maxmag));
line([leftx leftx+maxmag*sc]+tlen/20,maxmag*[1 1]*sc-yoff/2,'linestyle','-');
text(leftx+maxmag*sc+tlen/20,maxmag*sc-yoff/2,[sprintf(' %g ',maxmag) units]);

% labels
if ~isempty(labels)
  for i=1:M
    text(t(1),-yoff/2-(i-1)*yoff,labels(i,:),'horizontal','right');
  end;
end

set(gca,'Ytick',[]);
set(gca,'box','off');
% Warning: axes AspectRatio has been superseded by DataAspectRatio
% and PlotBoxAspectRatio and will not be supported in future releases.
%set(gca,'Aspect',[NaN 1]);
%set(gca,'DataAspectRatio',[NaN 1]);

%if (autoaxis), axis([min(tt)-max(wx) min(tt)+tlen+max(wx) -yax 0]);
%else axis([leftx leftx+tlen -yax 0]); end;
if (autoaxis), set(gca,'xlim',[min(tt)-max(wx) min(tt)+tlen+max(wx)]);
else set(gca,'xlim',[leftx leftx+tlen]); end;
% scale bar

%axis image %%% e



