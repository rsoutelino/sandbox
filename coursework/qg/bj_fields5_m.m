clear all ; close all ; clc

%**********************************

% **compute Velocity and N2 profiles
%   for BJ88 model**
%
%   VERSION FOR THE CRUISE AVERAGE 
%            45-51-60 


% modifieds to consider the current
% core bounded by the .15m/s isotach

%           October 2003

%**********************************

% load files

%cruz=[45 51 60];
%  cruz=[45 51];
%  
%  %n2=zeros(3,997);
%  n2=zeros(2,997);
%  v=n2;
%  
%  for i=1:2,
%  Cr=num2str(cruz(i));
%  ldfile=['load rperf',Cr];
%  eval(ldfile)
%  
%  n2(i,:)=n2r;
%  v(i,:)=vr';
%  end
%  
%  clear n2r vr
%  
%  vm=mean(v); stdv=std(v);
%  n2m=mean(n2); stdn2=std(n2);
%  pm=pr;
%  
%  vp=vm+stdv; vn=vm-stdv;
%  n2p=abs(n2m+stdn2); n2n=abs(n2m-stdn2);
%  
%  OP=['m'];

load N2.mat;

n2=N2;
p=p';

% *************************

%        N2 smoothing

% *************************

% smoothing begins
% polynomial fit on log(N2)
% offset is asymptote value of the N2 profile
% at the bottom



pix=[5:5:995];
pfit=binav(pix,p,p);
n2fit=binav(pix,p,n2);

%pfit(1)=10;
z=-pfit/1000;
N2=n2fit*86400^2;
bbar=-cumsum(N2)*0.005;

% ord=polynomial order 

ord=16;

offset=4e5; % for n2m 


phi=log(bbar+offset);

%x=0*z+1;
y=0*z;

for i=1:ord
%  x=[x,z.^i];
  y=[y,i*z.^(i-1)];
end

[beta,res]=polyfit(z,phi,ord);
phiest=polyval(beta,z);

%plot(z,phi,z,phiest)

n2=exp(phiest).*(y*fliplr(beta)');
n2x=exp(phi).*(y*fliplr(beta)');
bbarest=-offset+exp(phiest);

n2=n2/86400^2;
N2=N2/86400^2;

% fix the very upper part of the profile

n2f=n2;

% for n2m
n2f=n2f(4:length(n2f)); n2f=[ 1.e-5 ; 6e-5; n2f];
pfit=pfit(4:length(pfit)); pfit=[0; 10; pfit]; 

%n2f=n2f(2:length(n2f));  %n2f=[ 3.e-5 ; 6e-5; n2f];
%pfit=pfit(3:length(pfit)); pfit=[0; pfit]; 

figure
plot(n2f,-pfit,'r.','linewidth',2); hold on
plot(n2m,-pm,'b-')
hold off


% interpolate to standard depths

dzz=5;
zzi=-(5:dzz:895);

pn2=csape(-pfit,n2f,'variational');
n2=fnval(pn2,zzi);
mn2=mean(n2);

figure
plot(n2,zzi,'r.','linewidth',2); hold on
plot(n2m,-pm,'b-')
hold off


% differentiate function

pn=csape(-pfit,sqrt(n2f),'variational');
dn=fnder(pn,1);  
nz=fnval(dn,zzi);

dp=pfit(3)-pfit(2);
Nz=-gradient(sqrt(N2),dp);

%fix again

nz=weim(5,'hamm',nz);
n2=(weim(5,'hamm',sqrt(n2))).^2;


% ******figure 2********

tit=['Cruise Average N^2']
%prfile=['print -depsc modeln2_',OP,'.eps']

figure
subplot(121)
plot(1e4*n2,zzi,'r','linewidth',2)
axis([0 3.5 -900 0])
hold on
plot(1e4*N2,1000*z,'b')
hold off
title(tit)
xlabel('[10^{-4} rad^2 s^{-2}]')
ylabel('depth in meters')

subplot(122)
plot(1e4*nz,zzi,'m','linewidth',2)
axis([-4 4 -900 0])
hold on
plot(1e4*Nz,1000*z,'b')
plot([0 0],[0 -900],'k')
hold off
title('N_z')
xlabel('[10^{-4} rad (m s)^{-1}]')
ylabel('depth in meters')

drawnow
%eval(prfile)


%*****************************

% Velocity Profile

%******************************

global f;
global xi0;
global po;
global n2f2;

lV=length(zzi);

lat=-22;
f0= 14.585e-5*sin(lat*pi/180);
f2=f0*f0;
po=-(zzi(2:lV)+zzi(1:lV-1))/2;

n2f2=n2/f2;


rd=[25 16 8 6.5 5 4.5 3.5];

eigm=zeros(length(po),7);

for i=1:7,
r(i)=fzero('vmode2',rd(i))
eigm(:,i)=f';

%figure
%plot(eigm(:,i),-po,'b',[0 0],[-900 0],'k')
%pause

end



eigm=[ones(length(po),1) eigm];

F0=eigm(:,1);
F1=eigm(:,2);
F2=eigm(:,3);
F3=eigm(:,4);
F4=eigm(:,5);
F5=eigm(:,6);

zer = [0 0; 0 -900];

% ********figure 3**********

tit3=['Cruise Average Dynamical Modes']
prfile=['print -depsc modelmod_',OP,'.eps']

figure
plot(F0,-po,'r-.', F1,-po,'b-', F2,-po,'g--','linewidth',2)
axis([-5 6 -900 0])
axis('square')
hold on
plot(F3,-po,'m',F4,-po,'y--',F5,-po,'c-.','linewidth',2 )
legend('Zeroth','First', 'Second','Third','Fourth','Fifth',0)
plot(zer(:,1),zer(:,2),'k-')
hold off

title(tit3)
xlabel('Pressure Eigenmode Amplitudes')
ylabel('Depth in Meters')


orient tall
%eval(prfile)

% compute modal amplitudes

V=vm(11:5:901);
VV=(V(2:lV)+V(1:lV-1))/2;
dz=po(2)-po(1);
H=po(lV-1)-po(1);

for i=1:8,
  a(i)=1/H*sum(VV'.*eigm(:,i)*dz);
end

disp('the modal amplitudes [in m/s] are:')
[(0:7)'  a']

% compute modal rms

rms=zeros(5,1);

% for 1st baroclinic only

rms(1)= sqrt(mean((VV' - a(2)*eigm(:,2)).^2))/sqrt(mean(VV.^2));

% for 2nd baroclinic only

rms(2)= sqrt(mean((VV' - a(3)*eigm(:,3)).^2))/sqrt(mean(VV.^2));

% for 1st 2 modes

rms(3)= sqrt(mean((VV' -a(1)-a(2)*eigm(:,2)).^2))/sqrt(mean(VV.^2));

% for 1st 3 modes

rms(4)= sqrt(mean((VV' - a(1)-a(2)*eigm(:,2)-a(3)*eigm(:,3)).^2))/sqrt(mean(VV.^2));

% for 1st 4 modes

rms(5)= sqrt(mean((VV' - a(1)-a(2)*eigm(:,2)-a(3)*eigm(:,3)-a(4)*eigm(:,4)).^2))/sqrt(mean(VV.^2));

% reconstitute V profile with nmodes

nmode=input('How many modes in truncation: ');
nmode=4;
Vt=eigm(:,1:nmode)*a(1:nmode)';


% fix the very upper part of the vel profile


%Vf=[VV(1:3)'; Vt(10:end)];
%pp=[po(1:3)';po(10:end)'];

%pv=csape(-pp,Vf,'variational');
%Vt=fnval(pv,-po);

figure
plot(Vt,-po,'b-',VV,-po,'r.')

Vt=weim(3,'hamm',Vt)';
Vtz=-gradient(Vt,dz);

VVz=-gradient(VV,dz);


% ****figure 4******

TXT=[num2str(nmode),'-mode approx.']
tit2=['Cruise Average Meridional Velocity']
prfile=['print -depsc2 /ipa1/ilson/texstuff/projs/petro/rel4/modelV_',OP,'.eps']

figure
subplot(121)
plot(Vt,-po,'r','linewidth',2);hold on
plot(VV,-po,'b'); hold off
axis([-0.4 0.2 -900 0])
hold on
plot([0 0],[0 -900],'k')
hold off
title(tit2)
xlabel('[m/s]')
ylabel('depth in meters')
text(-0.35,-800,TXT)

subplot(122)
plot(1e3*Vtz,-po,'m','linewidth',2);hold on;
plot(1e3*VVz,-po,'c');hold off
axis([-3 2 -900 0])
hold on
plot([0 0],[0 -900],'k')
hold off
title('dV/dz')
xlabel('[1e3/s]')
ylabel('depth in meters')

% compute the potential vorticity zonal gradient profile Qx

n2m=(n2(2:lV)+n2(1:lV-1))/2;
nzm=(nz(2:lV)+nz(1:lV-1))/2;

Qx=-gradient(f0*f0*Vtz'./n2m,dz);
Qx(length(Qx))=Qx(length(Qx)-2);
Qx(length(Qx)-1)=Qx(length(Qx)-2);
%Qx=weim(5,'hamm',Qx);

figure

subplot(131)
plot(1e4*n2m,-po,'k','linewidth',2)
axis([0 3.5 -900 0])
%title(tit)
title('N^2(z)')
xlabel('[10^{-4} rad^2 s^{-2}]')
ylabel('depth in meters')
ylabel('profundidade [m]')

subplot(132)
plot(Vt,-po,'k','linewidth',2)
axis([-0.5 0.25  -900 0])
%plot(vs,zi,'b','linewidth',2)
hold on
plot([0 0],[0 -900],'k')
hold off
%title('Meridional Velocity')
title('Velocidade Meridional')
xlabel('[m/s]')
%ylabel('depth in meters')

subplot(133)
plot(1e9*Qx,-po,'k','linewidth',2)
axis([-0.5 1.2 -900 0])
hold on
plot([0 0],[0 -900],'k')
hold off
title('dq/dx')
xlabel('[10^{9}(m s)^{-1}]')
%ylabel('depth in meters')

%eval(prfile)
%print -djpeg90 /ipa1/ilson/bjfields.jpg

N2=n2m';Nz=nzm';V=Vt;Vz=Vtz;zi=-po';N=sqrt(N2);
Qx=Qx';
svfile=['save perf',OP,' N2 Nz N V Vz zi Qx'];
eval(svfile);
