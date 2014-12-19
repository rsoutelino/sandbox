%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INTERPOLA CONTORNOS HYCOM PARA ROMS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DADOS DO ROMS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=40;   %numero de camadas
theta_s=7.;
theta_b=0.;
hc=5.;

nc=netcdf('mesoscale_grd.nc');
lon_r=nc{'lon_rho'}(:);
lat_r=nc{'lat_rho'}(:);
lon_u=nc{'lon_u'}(:);
lat_u=nc{'lat_u'}(:);
lon_v=nc{'lon_v'}(:);
lat_v=nc{'lat_v'}(:);
h_roms=nc{'h'}(:);
close(nc)

im=length(lon_r(1,:));
jm=length(lon_r(:,1));
kb=N;

%%%%%%%%%%% CALCULO DE PROFUNDIDADES %%%%%%%%
zroms_r=squeeze(zlevs(h_roms,0.*h_roms,theta_s,theta_b,hc,N,'r'));
zroms_v=rho2v_3d(zroms_r);
zroms_u=rho2u_3d(zroms_r);

profr_south=squeeze(zroms_r(:,1,:));
profr_north=squeeze(zroms_r(:,end,:));
lonr_south =lon_r(1,:);
lonr_north =lon_r(end,:);
latr_south  =lat_r(1,:);
latr_north  =lat_r(end,:);

profr_east =squeeze(zroms_r(:,:,end));
profr_west =squeeze(zroms_r(:,:,1));
lonr_east  =lon_r(:,end);
latr_east  =lat_r(:,end);
lonr_west  =lon_r(:,1);
latr_west  =lat_r(:,1);

profu_south=squeeze(zroms_u(:,1,:));
profu_north=squeeze(zroms_u(:,end,:));
lonu_south =lon_u(1,:);
lonu_north =lon_u(end,:);
latu_south  =lat_u(1,:);
latu_north  =lat_u(end,:);

profu_east =squeeze(zroms_u(:,:,end));
profu_west =squeeze(zroms_u(:,:,1));
lonu_east  =lon_u(:,end);
latu_east  =lat_u(:,end);
lonu_west  =lon_u(:,1);
latu_west  =lat_u(:,1);

profv_south=squeeze(zroms_v(:,1,:));
profv_north=squeeze(zroms_v(:,end,:));
lonv_south =lon_v(1,:);
lonv_north =lon_v(end,:);
latv_south  =lat_v(1,:);
latv_north  =lat_v(end,:);

profv_east =squeeze(zroms_v(:,:,end));
profv_west =squeeze(zroms_v(:,:,1));
lonv_east  =lon_v(:,end);
latv_east  =lat_v(:,end);
lonv_west  =lon_v(:,1);
latv_west  =lat_v(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ABRINDO ARQUIVO BRY ROMS
nc=netcdf('FM_BC_bry.nc','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DADOS DO HYCOM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hycom_dir='saidas_2010/';
ano='2010'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DADOS DA GRADE DO HYCOM %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IDM=601; %numeros de pontod idm
JDM=733; %numeros de pontos jdm
KDM=21;  %numeros total de camadas
grid_fid=fopen('regional.grid.a','r'); %endereco do regional grid
IJDM=IDM*JDM;

SSH_INDEX=1;
UBARO_INDEX=6;
VBARO_INDEX=7;
UVEL_INDEX=8;
VVEL_INDEX=9;
LT_INDEX=10;
TEMP_INDEX=11;
SALT_INDEX=12;

NUM3DF=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npad=4096-mod(IJDM,4096);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%GRADE DO HYCOM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read plon and plat from regional grid file

[plon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
fseek(grid_fid,4*(npad+IJDM),-1);
[plat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

plon=reshape(plon,IDM,JDM);
plon=plon(:,1);
plat=reshape(plat,IDM,JDM);
plat=plat(1,:);

% read ulon and ulat from regional grid file
fseek(grid_fid,4*4*(npad+IJDM),-1);
[ulon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
fseek(grid_fid,5*4*(npad+IJDM),-1);
[ulat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

 ulon=reshape(ulon,IDM,JDM);
 ulon=ulon(:,1);
 ulat=reshape(ulat,IDM,JDM);
 ulat=ulat(1,:);
% 
% %%read vlon and vlat from regional grid file
fseek(grid_fid,6*4*(npad+IJDM),-1);
[vlon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
fseek(grid_fid,7*4*(npad+IJDM),-1);
[vlat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

 vlon=reshape(vlon,IDM,JDM);
 vlon=vlon(:,1);
 vlat=reshape(vlat,IDM,JDM);
 vlat=vlat(1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Lendo batimetria %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depth_fid=fopen('depth_ATLe0.08_03.a','r');

[batim,count]=fread(depth_fid,IJDM,'float32','ieee-be');
fseek(depth_fid,4*(npad+IJDM),-1);
y=find(batim>1e10);
batim(y)=nan;
batim=reshape(batim,IDM,JDM)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% COMECANDO LOOP NO TEMPO %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn=0;
%LOOP DOS DIAS
for n=1:10%:365   %defino os dias que serao usados
    nn=nn+1
    
    % monta nome do arquivo para cada tempo
    if n<10,            file='_00'; end
    if n>=10 & n<100,   file='_0'; end
    if n>=100 & n<1000, file='_'; end
    file1=['../HYCOM_CH_7d/saidas_2010/archv.',num2str(ano),file,num2str(n),'_00.a'];
    % aux=mod(n,30);if aux==0, file1, end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ABRINDO ARQUIVOS DO HYCOM
    
    archv_fid=fopen([hycom_dir,file1]);
    
    %%%%%% VARIAVEIS 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%%% SSH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fseek(archv_fid,SSH_INDEX*4*(npad+IJDM),-1);
    [ssh,count]=fread(archv_fid,IJDM,'float32','ieee-be');
    y=find(ssh>1e10);
    ssh(y)=nan;
    ssh=reshape(ssh,IDM,JDM)'.*1/(0.001*10*9.806);
    ssh=ssh./100;   %transformo p/ metros
    %%%%%% U BAROTROPICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fseek(archv_fid,UBARO_INDEX*4*(npad+IJDM),-1);
    [ubaro,count]=fread(archv_fid,IJDM,'float32','ieee-be');
    y=find(ubaro>1e10);
    ubaro(y)=nan;
    ubaro=reshape(ubaro,IDM,JDM)';
    %%%% V BAROTROPICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fseek(archv_fid,VBARO_INDEX*4*(npad+IJDM),-1);
    [vbaro,count]=fread(archv_fid,IJDM,'float32','ieee-be');
    y=find(vbaro>1e10);
    vbaro(y)=nan;
    vbaro=reshape(vbaro,IDM,JDM)';clear y;
    
    %%%%%%%% INTERPOLANDO HORIZONTALMENTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    zeta_south =interp2(plon,plat,ssh,lonr_south,latr_south);
    zeta_north =interp2(plon,plat,ssh,lonr_north,latr_north);
    zeta_east  =interp2(plon,plat,ssh,lonr_east,latr_east);
    zeta_west  =interp2(plon,plat,ssh,lonr_west,latr_west);
    ubar_south =interp2(ulon,ulat,ubaro,lonu_south,latu_south);
    ubar_north =interp2(ulon,ulat,ubaro,lonu_north,latu_north);
    ubar_east  =interp2(ulon,ulat,ubaro,lonu_east,latu_east);
    ubar_west  =interp2(ulon,ulat,ubaro,lonu_west,latu_west);
    vbar_south =interp2(vlon,vlat,vbaro,lonv_south,latv_south);
    vbar_north =interp2(vlon,vlat,vbaro,lonv_north,latv_north);
    vbar_east  =interp2(vlon,vlat,vbaro,lonv_east,latv_east);
    vbar_west  =interp2(vlon,vlat,vbaro,lonv_west,latv_west);
    
    %%%%%%%% TIRANDO TERRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    cc = find(isnan(zeta_west)==1); zeta_west(cc)=-1.e10;
    cc = find(isnan(zeta_east)==1); zeta_east(cc)=-1.e10;
    cc = find(isnan(zeta_south)==1); zeta_south(cc)=-1.e10;
    cc = find(isnan(zeta_north)==1); zeta_north(cc)=-1.e10;
    
    cc = find(isnan(ubar_west)==1); ubar_west(cc)=-1.e10;
    cc = find(isnan(ubar_east)==1); ubar_east(cc)=-1.e10;
    cc = find(isnan(ubar_south)==1); ubar_south(cc)=-1.e10;
    cc = find(isnan(ubar_north)==1); ubar_north(cc)=-1.e10;
    
    cc = find(isnan(vbar_west)==1); vbar_west(cc)=-1.e10;
    cc = find(isnan(vbar_east)==1); vbar_east(cc)=-1.e10;
    cc = find(isnan(vbar_south)==1); vbar_south(cc)=-1.e10;
    cc = find(isnan(vbar_north)==1); vbar_north(cc)=-1.e10;
    
    %%%%%%%%%%%%%% VARIAVEIS 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for k=0:KDM-1
        %%%%%%%%%% TEMPERATURA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fseek(archv_fid,(k*NUM3DF+TEMP_INDEX)*4*(npad+IJDM),-1);
        [a,count]=fread(archv_fid,IJDM,'float32','ieee-be');
        y=find(a>1e10);
        a(y)=nan;
        temp=reshape(a,IDM,JDM)'; clear a;
        clear y
        %%%%%%%%% ESPESSURA DAS CAMADAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fseek(archv_fid,(k*NUM3DF+LT_INDEX)*4*(npad+IJDM),-1);
        [a,count]=fread(archv_fid,IJDM,'float32','ieee-be');
        y=find(a>1e10);
        a(y)=nan;
        lthk=reshape(a,IDM,JDM)'./9806;clear a;
        clear y
        %%%%%%%%%% SALINIDADE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fseek(archv_fid,(k*NUM3DF+SALT_INDEX)*4*(npad+IJDM),-1);
        [a,count]=fread(archv_fid,IJDM,'float32','ieee-be');
        y=find(a>1e10);
        a(y)=nan;
        sal=reshape(a,IDM,JDM)'; clear a;
        clear y
        %%%%%%%%%%% U BAROCLINICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fseek(archv_fid,(k*NUM3DF+UVEL_INDEX)*4*(npad+IJDM),-1);
        [a,count]=fread(archv_fid,IJDM,'float32','ieee-be');
        y=find(a>1e10);
        a(y)=nan;
        u=reshape(a,IDM,JDM)'; clear a;
        %%%%%%%%%%% V BAROCLINICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fseek(archv_fid,(k*NUM3DF+VVEL_INDEX)*4*(npad+IJDM),-1);
        [a,count]=fread(archv_fid,IJDM,'float32','ieee-be');
        y=find(a>1e10);
        a(y)=nan;
        v=reshape(a,IDM,JDM)';clear a;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% INTERPOLANDO HORIZONTALMENTE
        
        temp_south(k+1,:) =interp2(plon,plat,temp,lonr_south,latr_south);
        sal_south(k+1,:)  =interp2(plon,plat,sal,lonr_south,latr_south);
        lthk_south(k+1,:) =interp2(plon,plat,lthk,lonr_south,latr_south);
        u_south(k+1,:)    =interp2(ulon,ulat,u,lonu_south,latu_south);
        v_south(k+1,:)    =interp2(vlon,vlat,v,lonv_south,latv_south);
        temp_north(k+1,:) =interp2(plon,plat,temp,lonr_north,latr_north);
        sal_north(k+1,:)  =interp2(plon,plat,sal,lonr_north,latr_north);
        lthk_north(k+1,:) =interp2(plon,plat,lthk,lonr_north,latr_north);
        u_north(k+1,:)    =interp2(ulon,ulat,u,lonu_north,latu_north);
        v_north(k+1,:)    =interp2(vlon,vlat,v,lonv_north,latv_north);
        temp_east(k+1,:)  =interp2(plon,plat,temp,lonr_east,latr_east);
        sal_east(k+1,:)   =interp2(plon,plat,sal,lonr_east,latr_east);
        lthk_east(k+1,:)  =interp2(plon,plat,lthk,lonr_east,latr_east);
        u_east(k+1,:)     =interp2(ulon,ulat,u,lonu_east,latu_east);
        v_east(k+1,:)     =interp2(vlon,vlat,v,lonv_east,latv_east);
        temp_west(k+1,:)  =interp2(plon,plat,temp,lonr_west,latr_west);
        sal_west(k+1,:)   =interp2(plon,plat,sal,lonr_west,latr_west);
        lthk_west(k+1,:)  =interp2(plon,plat,lthk,lonr_west,latr_west);
        u_west(k+1,:)     =interp2(ulon,ulat,u,lonu_west,latu_west);
        v_west(k+1,:)     =interp2(vlon,vlat,v,lonv_west,latv_west);
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% INTERPOLANDO VERTICALMENTE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% NORTH E SOUTH
    temp_south2=nan(kb,length(temp_south(1,:)));
    sal_south2 =nan(kb,length(temp_south(1,:)));
    u_south2 =nan(kb,length(u_south(1,:)));
    v_south2 =nan(kb,length(v_south(1,:)));
    temp_north2=nan(kb,length(temp_north(1,:)));
    sal_north2 =nan(kb,length(temp_north(1,:)));
    u_north2 =nan(kb,length(u_north(1,:)));
    v_north2 =nan(kb,length(v_north(1,:)));
    
    for i=1:im    % é im pois a grade está invertida
        
        if (isnan(temp_south(1,i))~=1),
            espessura1=lthk_south(:,i);
            temp_south1=interpz(espessura1,temp_south(:,i),-profr_south(:,i),KDM,1,1/(2.0001));
            sal_south1 =interpz(espessura1,sal_south(:,i),-profr_south(:,i),KDM,1,1/(2.0001));
            
            k=1;
            while (sum(espessura1)<-profr_south(k,i))
                temp_south1(1:k)=temp_south1(k+1); sal_south1(1:k)=sal_south1(k+1);
                k=k+1;
            end
            
            temp_south2(:,i)=temp_south1;
            sal_south2(:,i)=sal_south1;
            
        end
        
        if (isnan(temp_north(1,i))~=1),
            espessura1=lthk_north(:,i);
            temp_north1=interpz(espessura1,temp_north(:,i),-profr_north(:,i),KDM,1,1/(2.0001));
            sal_north1 =interpz(espessura1,sal_north(:,i),-profr_north(:,i),KDM,1,1/(2.0001));
            
            k=1;
            while (sum(espessura1)<-profr_north(k,i))
                temp_north1(1:k)=temp_north1(k+1); sal_north1(1:k)=sal_north1(k+1);
                k=k+1;
            end
            
             temp_north2(:,i)=temp_north1;
             sal_north2(:,i)=sal_north1;
             
        end        
         
    end
    
    for i=1:im-1
        
        if (isnan(v_south(1,i))~=1),
            espessura1=lthk_south(:,i);
            v_south1 =interpz(espessura1,v_south(:,i),-profv_south(:,i),KDM,1,1/(2.0001));
            
            k=1;
            while (sum(espessura1)<-profv_south(k,i))
                v_south1(1:k)=v_south1(k+1);
                k=k+1;
            end
            
            v_south2(:,i)=v_south1;
            
        end
        
        if (isnan(v_north(1,i))~=1),
            espessura1=lthk_north(:,i);
            v_north1 =interpz(espessura1,v_north(:,i),-profv_north(:,i),KDM,1,1/(2.0001));
            
            k=1;
            while (sum(espessura1)<-profv_north(k,i))
                v_north1(1:k)=v_north1(k+1);
                k=k+1;
            end
            
            v_north2(:,i)=v_north1;
        
        end        
        
    end
    
    for i=1:im-1
        
        if (isnan(u_south(1,i))~=1),
            espessura1=lthk_south(:,i);
            u_south1 =interpz(espessura1,u_south(:,i),-profu_south(:,i),KDM,1,1/(2.0001));
            
            k=1;
            while (sum(espessura1)<-profu_south(k,i))
                u_south1(1:k)=u_south1(k+1);
                k=k+1;
            end
            
            u_south2(:,i)=u_south1;
            
        end
        
        if (isnan(u_north(1,i))~=1),
            espessura1=lthk_north(:,i);
            u_north1 =interpz(espessura1,u_north(:,i),-profu_north(:,i),KDM,1,1/(2.0001));
            
            k=1;
            while (sum(espessura1)<-profu_north(k,i))
                u_north1(1:k)=u_north1(k+1);
                k=k+1;
            end
            
            u_north2(:,i)=u_north1;
            
        end   
                
    end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EAST E WEST
    temp_east2=nan(kb,length(temp_east(1,:)));
    sal_east2=nan(kb,length(temp_east(1,:)));
    u_east2=nan(kb,length(u_east(1,:)));
    v_east2=nan(kb,length(v_east(1,:)));
    temp_west2=nan(kb,length(temp_west(1,:)));
    sal_west2=nan(kb,length(temp_west(1,:)));
    u_west2=nan(kb,length(u_west(1,:)));
    v_west2=nan(kb,length(v_west(1,:)));
    
    
    for i=1:jm  %jm pois a grade está invertida
        
        if (isnan(temp_east(1,i))~=1),
            espessura1=lthk_east(:,i);
            temp_east1=interpz(espessura1,temp_east(:,i),-profr_east(:,i),KDM,1,1/(2.0001));
            sal_east1 =interpz(espessura1,sal_east(:,i),-profr_east(:,i),KDM,1,1/(2.0001));
                        
            k=1;
            while (sum(espessura1)<-profr_east(k,i))
                temp_east1(1:k)=temp_east1(k+1); sal_east1(1:k)=sal_east1(k+1);
                k=k+1;
            end
            
            temp_east2(:,i)=temp_east1;
            sal_east2(:,i)=sal_east1;
            
        end
        
            
        if (isnan(temp_west(1,i))~=1),
            espessura1=lthk_west(:,i);
            temp_west1=interpz(espessura1,temp_west(:,i),-profr_west(:,i),KDM,1,1/(2.0001));
            sal_west1 =interpz(espessura1,sal_west(:,i),-profr_west(:,i),KDM,1,1/(2.0001));
                        
            k=1;
            while (sum(espessura1)<-profr_west(k,i))
                temp_west1(1:k)=temp_west1(k+1); sal_west1(1:k)=sal_west1(k+1);
                k=k+1;
            end
            
            temp_west2(:,i)=temp_west1;
            sal_west2(:,i)=sal_west1;
        end
        
    end
    
    
    
    for i=1:jm  %jm pois a grade está invertida
        
        if (isnan(u_east(1,i))~=1),
            espessura1=lthk_east(:,i);
            u_east1 =interpz(espessura1,u_east(:,i),-profr_east(:,i),KDM,1,1/(2.0001));
                        
            k=1;
            while (sum(espessura1)<-profr_east(k,i))
                u_east1(1:k)=u_east1(k+1); 
                k=k+1;
            end
            
            u_east2(:,i)=u_east1;
            
        end
        
    end
    
    for i=1:jm-1  %jm pois a grade está invertida
        
        if (isnan(v_east(1,i))~=1),
            espessura1=lthk_east(:,i);
            v_east1 =interpz(espessura1,v_east(:,i),-profr_east(:,i),KDM,1,1/(2.0001));
                        
            k=1;
            while (sum(espessura1)<-profr_east(k,i))
                v_east1(1:k)=v_east1(k+1); 
                k=k+1;
            end
            
            v_east2(:,i)=v_east1;
            
        end
        
    end   
     
        
    for i=1:jm
        
        if (isnan(u_west(1,i))~=1),
            espessura1=lthk_west(:,i);
            u_west1 =interpz(espessura1,u_west(:,i),-profu_west(:,i),KDM,1,1/(2.0001));
                        
            k=1;
            while (sum(espessura1)<-profu_west(k,i))
                u_west1(1:k)=u_west1(k+1);
                k=k+1;
            end
            
           u_west2(:,i)=u_west1;
        end
        
    end
    
    for i=1:jm-1
        
        if (isnan(v_west(1,i))~=1),
            espessura1=lthk_west(:,i);
            v_west1 =interpz(espessura1,v_west(:,i),-profv_west(:,i),KDM,1,1/(2.0001));
                        
            k=1;
            while (sum(espessura1)<-profv_west(k,i))
                v_west1(1:k)=v_west1(k+1);
                k=k+1;
            end
            
            v_west2(:,i)=v_west1;
        end
        
    end   
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cc = find(isnan(temp_west2)==1); temp_west2(cc)=-1.e10;
    cc = find(isnan(temp_east2)==1); temp_east2(cc)=-1.e10;
    cc = find(isnan(temp_south2)==1); temp_south2(cc)=-1.e10;
    cc = find(isnan(temp_north2)==1); temp_north2(cc)=-1.e10;
    
    cc = find(isnan(sal_west2)==1); sal_west2(cc)=-1.e10;
    cc = find(isnan(sal_east2)==1); sal_east2(cc)=-1.e10;
    cc = find(isnan(sal_south2)==1); sal_south2(cc)=-1.e10;
    cc = find(isnan(sal_north2)==1); sal_north2(cc)=-1.e10;
    
    cc = find(isnan(u_west2)==1); u_west2(cc)=-1.e10;
    cc = find(isnan(u_east2)==1); u_east2(cc)=-1.e10;
    cc = find(isnan(u_south2)==1); u_south2(cc)=-1.e10;
    cc = find(isnan(u_north2)==1); u_north2(cc)=-1.e10;
    
    cc = find(isnan(v_west2)==1); v_west2(cc)=-1.e10;
    cc = find(isnan(v_east2)==1); v_east2(cc)=-1.e10;
    cc = find(isnan(v_south2)==1); v_south2(cc)=-1.e10;
    cc = find(isnan(v_north2)==1); v_north2(cc)=-1.e10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SALVANDO NO ARQUIVO DE BRY DO ROMS
    
    nc{'temp_north'}(n,:,:)=temp_north2;
    nc{'temp_south'}(n,:,:)=temp_south2;
    nc{'temp_east'}(n,:,:) =temp_east2;
    nc{'temp_west'}(n,:,:) =temp_west2;
    
    nc{'salt_north'}(n,:,:)=sal_north2;
    nc{'salt_south'}(n,:,:)=sal_south2;
    nc{'salt_east'}(n,:,:) =sal_east2;
    nc{'salt_west'}(n,:,:) =sal_west2;
    
    nc{'u_north'}(n,:,:)=u_north2;
    nc{'u_south'}(n,:,:)=u_south2;
    nc{'u_east'}(n,:,:) =u_east2;
    nc{'u_west'}(n,:,:) =u_west2;
    
    nc{'v_north'}(n,:,:)=v_north2;
    nc{'v_south'}(n,:,:)=v_south2;
    nc{'v_east'}(n,:,:) =v_east2;
    nc{'v_west'}(n,:,:) =v_west2;
    
    nc{'zeta_north'}(n,:)=zeta_north;
    nc{'zeta_south'}(n,:)=zeta_south;
    nc{'zeta_east'}(n,:) =zeta_east;
    nc{'zeta_west'}(n,:) =zeta_west;
    
    nc{'ubar_north'}(n,:)=ubar_north;
    nc{'ubar_south'}(n,:)=ubar_south;
    nc{'ubar_east'}(n,:) =ubar_east;
    nc{'ubar_west'}(n,:) =ubar_west;
    
    nc{'vbar_north'}(n,:)=vbar_north;
    nc{'vbar_south'}(n,:)=vbar_south;
    nc{'vbar_east'}(n,:) =vbar_east;
    nc{'vbar_west'}(n,:) =vbar_west;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FECHANDO ARQUIVO HYCOM
    fclose(archv_fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIM DO LOOP DO TEMPO
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FECHANDO ARQUIVO BRY HYCOM
close(nc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%