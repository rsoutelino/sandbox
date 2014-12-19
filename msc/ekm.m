load lat1.dat; load lon1.dat; load fsm1.dat
load wstrsx.dat; load wstrsy.dat
wxx=reshape(wstrsx',389,308,12); clear wstrsx
wyy=reshape(wstrsy',389,308,12); clear wstrsy
wxm=nanmean(wxx,3)'; wym=nanmean(wyy,3)';
ind=find(fsm1==0); wxm(ind)=NaN; wym(ind)=NaN;

load bif_dyn

wx=interp2(lon1,lat1,wxm,X,Y);
wy=interp2(lon1,lat1,wym,X,Y);

f=2*(2*pi/86400)*sin(Y*pi/180);
ind=find(f==0);f(ind)=NaN;

mex=wy./(f*1024);
mey=-wx./(f*1024);
de=pi*sqrt(0.01./abs(f)); dd=de/2;
%dy=ones(31,41)*111120;
%dx=cos(Y*pi/180)*111120;

ue=mex./dd; ve=mey./dd;
%ue=mex./20; ve=mey./20;
ind=find(Y>-5);  ue(ind)=NaN; ve(ind)=NaN; dd(ind)=NaN;
ind=find(Y<-27); ue(ind)=NaN; ve(ind)=NaN; dd(ind)=NaN;


in=~isnan(v(:,:,1));ind=cumsum(in,2); vm=NaN*ones(31,10);
for j=1:31
ind1=find(ind(j,:)>=1 & ind(j,:)<=5);
if isempty(ind1)
    vm(j,:)=NaN;
else
    vm(j,:)=squeeze(nanmean(v(j,ind1,:),2))';
end
end

Z=[0:100:900];contour(Y(:,1),Z',vm'*10,[-0.2:0.02:0.2],'-b'); axis ij
hold on; contour(Y(:,1),Z',vm'*10,[0 0],'-r')
axis([-30 -5 0 1000]);grid on;

vv=ve; %vv=v;
%vv(:,:,1)=v(:,:,1)+ve(:,:);
%vv(:,:,2)=v(:,:,2)+ve(:,:);
in=~isnan(vv(:,:,1));ind=cumsum(in,2); vm1=NaN*ones(31,1);
for j=1:31
ind1=find(ind(j,:)>=1 & ind(j,:)<=5);
if isempty(ind1)
    vm1(j,:)=NaN;
else
    vm1(j,:)=squeeze(nanmean(vv(j,ind1,:),2))';
end
end

vmm=vm*10; vmm(:,1)=vmm(:,1)+vm1;
figure(2)
Z=[0:100:900];contour(Y(:,1),Z',vmm',[-0.2:0.02:0.2],'-b'); axis ij
axis([-30 -5 0 1000]);grid on;hold on
contour(Y(:,1),Z',vmm',[0 0],'-r')
