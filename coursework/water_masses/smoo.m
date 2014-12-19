%***********************************************************
%
% Programa para tratamento de dados de CTD projeto ABROLHOS
%                    Campanha ABROLHOS 1
%
%***********************************************************


clear;close all;clc;

Se=33:.25:38;
Te=0:30;

[Sg,Tg]=meshgrid(Se,Te);
dens=sw_dens0(Sg,Tg)-1000;
xlabel('Salinidade [ups]')
     ylabel('Temperatura [graus C]')

figure(1)
     contour(Se,Te,dens,20:.2:40,'k')
hold on
drawnow

xx=1:12; %input(' Entre vetor com as estacoes que deseja filtrar.  Exemplo : [7550:75...]: ');
lx=length(xx);

for k=1:lx;
XX=num2str(xx(k));

% ---------------------------------------

file=['load dados/est',XX,'.dat'];
name=['est',XX];

eval(file);

ctd=eval(name);

% removendo spikes ----------------------

%  for i=1:size(ctd,1);
%  	if ctd(i,2) < 0;
%  	   ctd(i,:) = ctd(i-1,:);
%  	end
%  end	    

%  file2=['save -ascii ctd doei0',XX,'.dat'];
%  eval(file2);

end

 

p=ctd(1:end,1);
T=ctd(1:end,2);
S=ctd(1:end,3);

% if xx==7552
%    p=flipud(ctd(:,2));
%    T=flipud(ctd(:,3));
%    C=flipud(ctd(:,4));
%    SS=flipud(ctd(:,10));
% end

% converter condutividade para salinidade

%  R=1/sw_c3515*C;
%  S=sw_salt(R,T,p);


maxp=max(p); minp=min(p);
pmin=ceil(minp)+1;
pmax=floor(maxp)-1;

f=find(p>=pmin & p<=pmax); 
p=p(f);
T=T(f);
S=S(f);
%  SS=SS(f);
%  CC=C(f);
l1=length(p);

eval(['clear ',name]);


% fazer media em caixa (bin average)


binp=pmin:pmax;
lp=length(binp);

 for kk=1:lp-1
 	f=find(p>=binp(kk) & p<=binp(kk+1));
        Tbin(kk)=mean(T(f));%CCbin(kk)=mean(CC(f));
        Sbin(kk)=mean(S(f));%SSbin(kk)=mean(SS(f));
        pbin(kk)=.5*(binp(kk)+binp(kk+1));
end

% aplicar janela movel

if pmax>=500,
	Ts=weim(25,'hann',Tbin);  %CCs=weim(25,'hann',CCbin);
        Ss=weim(25,'hann',Sbin);  %SSs=weim(25,'hann',Sbin);
elseif pmax>=100 & pmax<500
	Ts=weim(11,'hann',Tbin); %CCs=weim(11,'hann',CCbin);
        Ss=weim(11,'hann',Sbin); %SSs=weim(11,'hann',Sbin);
else
	Ts=weim(5,'hann',Tbin); %CCs=weim(5,'hann',CCbin);
        Ss=weim(5,'hann',Sbin); %SSs=weim(5,'hann',Sbin);
end


t1=['Estacao ' ,XX,'--Temperatura com \Delta p = 1 dbar'];
t2=['Estacao ' ,XX,'--Salinidade com \Delta p = 1 dbar'];


% ----- correcao do desvio de -0.28 psu ----- %

%Ss=Ss-0.28;

% ------------------------------------------- %

figure(1)
     plot(Ss,Ts,'r.')
     hold on; %plot(Ss,Ts,'b.'); 
     title('Diagrama T-S espalhado')
drawnow

figure
plot(Ts,-pbin,'r','linewidth',3);
%plot(Tbin,-pbin,'r','linewidth',3);
hold on
plot(T,-p,'b--')
axis([0 30 -pbin(end) -1])
hold off
title(t1)
xlabel('Temperatura [graus C]')
ylabel('profundidade [metros]')
drawnow

figure
plot(Ss,-pbin,'g','linewidth',3);
%plot(Sbin,-pbin,'r','linewidth',3);
hold on
plot(S,-p,'b--')
axis([33 38 -pbin(end) -1])
hold off
title(t2)
xlabel('Salinidade [ups]')
ylabel('profundidade [metros]')
drawnow

%disp(['pause'])
%pause

     T=Ts'; S=Ss'; p=pbin'; C=CCs';


% calculo de T-Theta
TTH=sw_ptmp(S,T,p,0);

%calculo de sigma-T
sigT=sw_dens0(S,T)-1000;

%calculo de sigma-Theta
sigTH=sw_pden(S,T,p,0)-1000;

%calculo da velocidade do som
svel = sw_svel(S,T,p);
format long g
data=[p T C TTH S sigT sigTH svel]';


%  RR2=['fid =fopen(''doei',XX,'.dat'',''w'');'];

%  eval(RR2)

%         fprintf(fid,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n',data);


%  clear p T S C CC CCbin CCs binp pbin Tbin Sbin Ss Ts pmax pmin maxp minpl1 lp lx t1 t2

end

%figure(1)

saveas(gcf,'TS','fig');
