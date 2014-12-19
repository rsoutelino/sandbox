%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   WORLD OCEAN ATLAS
% Plota seções verticais de velocidade via metodo dinamico
% e salva matrizes de gpan para usar posteriormente nos 
%                    calculos de psi
%      Rafael Guarino Soutelino - Mestrado IOUSP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;clc;

%%% LIMITES OESTE-LESTE e SUL-NORTE DA AREA DE INTERESSE %%%

lonlim = [-41 -33];
latlim = [-20 -10];

%%% PERIODO DESEJADO %%%

disp(' ');
disp(['  ### Escolha do Periodo ###']);
disp(' ');
disp(['  - Para o campo anual: [anual]']);
disp(['  - Para o campo mensal: [jan] ou [fev] ou ...']);
disp(['  - Para media dos campos mensais: [jan,fev,mar,...]']);
disp(' ');
time = input(['  Periodo: '],'s');

time([1 end]) = [];

ft = find(time==',');

if isempty(ft) == 1
  Time = cellstr(time);
else
  for i = 1:length(ft)
    Time{i} = time(ft(i)-3:ft(i)-1);
  end
  Time{i+1} = time(ft(end)+1:ft(end)+3);
end


%%% LEITURA DA CLIMATOLOGIA DECLARADA %%%

for i = 1:length(Time)

  eval(['load  ../../../dados/WOA2001/',Time{i},'/Twoa_',Time{i},'_brazilcoast.dat;']);
  eval(['load  ../../../dados/WOA2001/',Time{i},'/Swoa_',Time{i},'_brazilcoast.dat;']);

  eval(['Twoa = Twoa_',Time{i},'_brazilcoast;']);
  eval(['Swoa = Swoa_',Time{i},'_brazilcoast;']);
  
  eval(['clear Twoa_',Time{i},'_brazilcoast;']);
  eval(['clear Swoa_',Time{i},'_brazilcoast;']);
  
  lat = Twoa(:,1)';
  lon = Twoa(:,2)';
  
  flat = find(lat >= min(latlim) & lat <= max(latlim));
  flon = find(lon >= min(lonlim) & lon <= max(lonlim));


  ii = 1;
  for jj = 1:length(flon)
    if isempty(find(flat == flon(jj))) == 0
      q(ii) = flon(jj);
      ii = ii+1;
    end
  end

  eval(['T',Time{i},' = Twoa(q,3:end);']);
  eval(['S',Time{i},' = Swoa(q,3:end);']);

end

%%% MEDIA DOS CAMPOS DESEJADOS %%%

strT = [];
strS = [];
for i = 1:length(Time)
  strT = [strT,'T',Time{i},'+'];
  strS = [strS,'S',Time{i},'+'];
end
strT(end) = [];
strS(end) = [];

eval(['Twoam = ((',strT,')/',num2str(i),')'';']);
eval(['Swoam = ((',strS,')/',num2str(i),')'';']);

latw = lat(q);
lonw = lon(q);

%%% MANIPULACAO DOS CAMPOS %%%

if size(time,2) == 5
  Pw = [0;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;...
        800;900;1000;1100;1200;1300;1400;1500;1750;2000;2500;3000;...
        3500;4000;4500;5000;5500];
else
  Pw = [0;10;20;30;50;75;100;125;150;200;250;300;400;500;600;700;...
        800;900;1000;1100;1200;1300;1400;1500];
end

pplot = 3000; % profmax
pi= % profundidade de interesse para plotar psi
pf = near(Pw,pplot,1);
z = Pw(1:pf);

clear Twoa Swoa Sanual Tanual f ft jj ii q strT strS


%%% MONTANDO MATRIZES DE LAT LON
% descobrindo tamanho da grade

f = find(diff(lonw)~=0);
li = f(2)-f(1);
lj = length(f)+1;
lon = reshape(lonw,li,lj);
lat = reshape(latw,li,lj);

gpan = [];

%começa looping para as radiais WOA
%  for k=1:li	

	for j=1:lj
           lo = lon(k,j); 
           la = lat(k,1);
           flo = find(lonw==lo);
           f = flo(1);
           T(:,j) = Twoam(1:pf,f);
	   S(:,j) = Swoam(1:pf,f);
 	   lons = lon(1,:);
           lats = lat(1,:);
	end

	[lz,lx]=size(T);
	zmax=max(z);
	[dst,ang]=sw_dist(lats,lons,'km'); % calculando distancia entre as estacoes
	x=[0 cumsum(dst)];
	xmax=max(x);

	clear ang f dst

        Ti=T; Si=S; xi=x;

   disp('Extrapolando campos')
   decT = 2; % taxa de decaimento do gradiente
   decS = 2;

   for i = 1:length(z)   
      f = find(isnan(Ti(i,:))==1); % descobrindo os nan da linha
      if isempty(f)==1; clc
      else
        ff = find(diff(f) > 1);
        if isempty(ff)==1  % teste para saber se ha montes submarinos

           % extrapolando o talude num perfil sem monte
           if f(end)+2 > lx % checando se ha alguem para calcular o gradiente
               for k=fliplr(f)
                 Tult = Ti(i,k+1); Sult = Si(i,k+1);
                 Ti(i,k) = Tult; Si(i,k) = Sult;
               end
           else
               dt = Ti(i,f(end)+2) - Ti(i,f(end)+1); % grad entre as 2 ult. estacoes   
               ds = Si(i,f(end)+2) - Si(i,f(end)+1);
               for k=fliplr(f)
                 dt = dt/decT; ds = ds/decS;
                 Tult = Ti(i,k+1); Sult = Si(i,k+1);                      
                 Tex = Tult - dt; Sex = Sult - ds;        
                 Ti(i,k) = Tex; Si(i,k) = Sex;              
               end
           end

        else

            % extrapolando o talude num perfil com monte
            f1 = f(1:ff); % NaNs do talude
            if isnan(Ti(i,ff+2)) == 1 % checando se ha alguem para calcular o gradiente
               for k=fliplr(f1)
                 Tult = Ti(i,k+1); Sult = Si(i,k+1);
                 Ti(i,k) = Tult; Si(i,k) = Sult;
               end
            else             
               dt = Ti(i,f1(end)+2) - Ti(i,f1(end)+1); % grad entre as 2 ult. estacoes   
               ds = Si(i,f1(end)+2) - Si(i,f1(end)+1);
               for k=fliplr(f1)
                 dt = dt/decT; ds = ds/decS;
                 Tult = Ti(i,k+1); Sult = Si(i,k+1);                      
                 Tex = Tult - dt; Sex = Sult - ds;        
                 Ti(i,k) = Tex; Si(i,k) = Sex;              
               end     
            end

            % extrapolando lado esquerdo do monte
            f2 = f(ff+1:end); % NaNs do monte
         
            for k=f2(1:ceil(length(f2)/2))
              Tult = Ti(i,k-1); Sult = Si(i,k-1);                             
              Ti(i,k) = Tult; Si(i,k) = Sult;              
            end   
            % extrapolando lado direito do monte
            for k=fliplr(f2(ceil(length(f2)/2):end))
              Tult = Ti(i,k+1); Sult = Si(i,k+1);                              
              Ti(i,k) = Tult; Si(i,k) = Sult;              
            end  
               
         end
      end
   end

   clear i k f ff f1 f2 ds dt 

	Di=sw_dens0(Si,Ti)-1000;

        %%% CALCULANDO ANOMALIA DO GEOPOTENCIAL

	nr = 1000;
	f0 = sw_f(-16); % latitude central da grade OEII
	% anomalia do geopotencial relativa a superficie
	agp = sw_gpan(Si,Ti,z);
	
	% anomalia do geopotencial relativa ao NR
	agp = agp - ones(size(z)) * agp(nr,:);
	
	gpan = [gpan]

      







