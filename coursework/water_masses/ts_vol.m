%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Diagrama TS Volumetrico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all;clc

camp = ('Lista 2');
%  r=input('Entre com a radial: Sul (s); Norte (n): ','s')
n_est=('est');
est=1:12;

load posicoes.mat


dist_1=sw_dist([lat(1) lat(6)],[lon(1) lon(6)],'Km');
dist_2=sw_dist([lat(1) lat(12)],[lon(1) lon(12)],'Km');
reg='Radial Norte';

area=dist_1*dist_2;

l_est=length(est);

nmax = zeros(l_est,1);

%compatibilizando as matrizes considerando a maior profundidade 

Tn=ones(3000,l_est)*NaN;
Sn=ones(3000,l_est)*NaN;
Pn=ones(3000,l_est)*NaN;

%file=['load ',num2str(lat),'.dat'];
%eval(file);
%file=['load ',num2str(lon),'.dat'];
%eval(file);

for k=1:l_est
file=['load dados/',num2str(n_est),num2str(est(k)),'.dat'];
eval(file); 
fname=[num2str(n_est),num2str(est(k))];
a=eval(fname);

% DADOS AMOSTRADOS
Sn(1:length(a),k)=a(:,3);
Tn(1:length(a),k)=a(:,2);
Pn(1:length(a),k)=a(:,1);

name=[num2str(n_est),num2str(est)];

end


dt=input(' Entre com o intervalo de variacao da Temperatura: ');
ds=input(' Entre com o intervalo de variacao da Salinidade: ');

[i,j]=size(Sn);
area_p=area/j;
maxc=j;




ii=zeros(30,1)*nan;
PP=zeros(3000,1);

total=0;


for l=1:maxc

i=find(isfinite(Sn(:,l)));
S=Sn(i,l);
T=Tn(i,l);
P=Pn(i,l);

prof(1,l)=P(end);


if (P(1,1) ~= 0)

  P(1,1)=0;

end



i=1;

soma=0;
soma_T=0;

n=0;
tmax=28;
smax=38;
smin=34;
tmin=0;

for ti=tmax:-2*dt:tmin
T_var=[ti ti-2*dt]
n=n+1;

it=find((T <= ti & T >= ti-2*dt));
SS=nan*S;

SS(it)=S(it);
soma_T=0;
k=0;

for si=smax:-2*ds:smin


%  Tm=mean(T_var)
%  Sm=mean(S_var)
  
k=k+1;
soma_pS=0;
is=find((SS < si & SS >= si-2*ds))


 if (isempty(is) == 0)
   T_var=[ti ti-2*dt]
   S_var=[si si-2*ds]
  if is(1) == 1
  a=[P(2) P(1)];
  dz=[P(2)-P(1)];
  %a=[P(is(2):is(length(is))) P(is(1):is(length(is))-1)];
  %dz=sum(P(is(2):is(length(is)))-P(is(1):is(length(is))-1));  
  else
  a=[P(is(1)-1:is(length(is))-1) P(is(1):is(length(is)))];
  dz=sum(P(is(1):is(length(is)))-P(is(1)-1:is(length(is))-1));
  end
  soma_pS=soma_pS+dz;
  soma=dz.*area*1e-3
%  Tm=mean(T_var);
%  Sm=mean(S_var);

%aux(n,k)=soma_pS

end

%Sm(n,k)=Sm;
aux(n,k)=soma_pS.*area_p*1e-3;
T(is)=nan;
S(is)=nan;
soma_T=soma_T+soma_pS;
end

soma_T=soma_T;
ST(n)=soma_T;

end

total=total+aux;

aux=aux.*0;
end

Prof_med=mean(prof);


figure(1)
orient landscape
t=(tmax:-dt:tmin)';
s=(smax:-ds:smin);
[Tg,Sg]=meshgrid(t,s);
sigmat=sw_dens0(Sg,Tg)-1000;
c=contour(Sg,Tg,sigmat,'k-');
clabel(c,'manual')
%set(gca,'plotboxaspectratio',[2 1 1])
set(gca,'YTick',[tmin:2*dt:tmax])
set(gca,'XTick',[smin:2*ds:smax])

grid
hold on


s=smax:-2*ds:smin;
t=(tmax:-2*dt:tmin)';
[i,j]=size(total);
tt=t*ones(1,j);
ss=ones(i,1)*s;
ii=find(total);
tot=total(ii);
ss=ss(ii);
tt=tt(ii);


%i_t=find(ss >= 37 & tt > 24);
i_t=find(ss >= 36 & tt > 20);
v_at=sum(tot(i_t));
%i_c=find((37.2 >= ss & ss >= 35.2) & (24 >= tt & tt > 12.));
i_c=find((36.4 >= ss & ss >= 35.2) & (20 >= tt & tt > 12.));
v_acas=sum(tot(i_c));
%i_a1=find((ss >= 34.2 & ss <= 35.2) & (12. >= tt & tt > 4));
%i_a2=find((ss >= 34.2 & ss <= 34.6) & (tt<=4));
i_a1=find((ss >= 34.2 & ss <= 35.2) & (12. >= tt & tt > 4));
i_a2=find((ss >= 34.2 & ss <= 34.6) & (tt<=4));
v_aia=sum([tot(i_a1) ;tot(i_a2)]);
i_ap=find((ss > 34.6 & ss < 35.2) & (tt <=4 & tt >1));
v_apan=sum(tot(i_ap));

Vtotal=v_at+v_acas+v_aia+v_apan;
v_at_p=v_at*100/Vtotal;
v_acas_p=v_acas*100/Vtotal;
v_aia_p=v_aia*100/Vtotal;
v_apan_p=v_apan*100/Vtotal;

Volume=[[v_at ;v_acas ;v_aia; v_apan] [v_at_p ;v_acas_p ;v_aia_p; v_apan_p]]  
%  eval(['save -ascii vol_tab_',r,'.dat  Volume'])

a=sprintf('%3.2f\n',tot);
a=str2num(a);

testo=text(ss-2*ds,tt-dt,num2str(a),'Color','b')
set(testo,'BackgroundColor',[.7 .9 .7])
%text(ss-ds,tt-dt,num2str(tot),'Color','r')
%text(ss(i_t)-(ds+ds/2),tt(i_t)-dt,num2str(tot(i_t)),'Color','r')
%text(ss(i_c)-(ds+ds/2),tt(i_c)-dt,num2str(tot(i_c)),'Color','m')
%text(ss(i_a)-(ds+ds/2),tt(i_a)-dt,num2str(tot(i_a)),'Color','b')
%text(ss(i_ap)-(ds+ds/2),tt(i_ap)-dt,num2str(tot(i_ap)),'Color','k')


str1(1)={['AT=',sprintf('%2.2f',v_at_p),'%']}
str1(2)={['ACAS=',sprintf('%2.2f',v_acas_p),'%']}
str1(3)={['AIA=',sprintf('%2.2f',v_aia_p),'%']}
str1(4)={['APAN=',sprintf('%2.2f',v_apan_p),'%']}
%str1(1)={['AT=',num2str(v_at_p), '%']}
%str1(2)={['ACAS=',num2str(v_acas_p), '% ']}
%str1(3)={['AIA=',num2str(v_aia_p), '%']}
%str1(4)={['APAN=',num2str(v_apan_p), '%']}
%uicontrol('style','text','Position',[440 150 130 60],'string',str1)
uicontrol('style','text','Position',[440 120 130 70],'string',str1,'FontSize',12,'FontWeight','demi')
title({['TS Volumerico Estatistico','FontSize',12,'FontWeight','demi')
xlabel('Salinidade','FontSize',12)
ylabel('Temperatura [^{\circ}{\rmC}]','FontSize',12)

%fn=['bio2_wb',r,'_TS_vol'];
fn=['sueli1'];
%  print( gcf, '-djpeg99 ', fn )
%  saveas(gcf,fn, 'fig');
%  eval(['print -depsc2 ',fn])



figure(2)
%  load ts_colormap
%  colormap(ts_colormap)
surfc(s,t,total)
shading interp
caxis([0 125]);
axis([34 38 0 28 0 200])
hold on
%mesh(s,t,total,'k')
view(20.5,54)
%set(gca,'XTick',[0 2 4 6 8 10 12 14 16 18 20])
%set(gca,'XTickLabel',[38.4 38 37.6 37.2 36.8 36.4 36 35.6 35.2 34.8 34.4 34 33.8])
%set(gca,'YTick',[0:4:30])
%set(gca,'YTickLabel',[28  24 20 16  12  8  4  0])
xlabel('Salinidade','FontSize',12)
ylabel('Temperatura [^{\circ}{\rmC}]','FontSize',12)
%title({['TS Volum�rico Estat�tico -',camp];[reg]},'FontSize',12,'FontWeight','demi')
zlabel('Volume','FontSize',12)
str1={['AT=',sprintf('%2.2f',v_at_p),'%']}
str2={['ACAS=',sprintf('%2.2f',v_acas_p),'%']}
str3={['AIA=',sprintf('%2.2f',v_aia_p),'%']}
str4={['APAN=',sprintf('%2.2f',v_apan_p),'%']}
text(35,24,200,str1,'Color',[0 0 0],'Fontweight','bold','FontSize',16)
text(35,24,160,str2,'Color',[0 0 0],'Fontweight','bold','FontSize',16)
text(35,24,120,str3,'Color',[0 0 0],'Fontweight','bold','FontSize',16)
text(35,24,80,str4,'Color',[0 0 0],'Fontweight','bold','FontSize',16)
fn=['sueli2'];
%  print( gcf, '-djpeg99 ', fn )
%  saveas(gcf,fn, 'fig');
%  eval(['print -depsc2 ',fn])
