clear all;close all;clc

for k = 32:92
	eval(['[lon,lat,u] = read_qscat_stress(''QS_STGRD3_20050',num2str(k),''',''asc_stress_Liu_U'');']);
	eval(['[lon,lat,v] = read_qscat_stress(''QS_STGRD3_20050',num2str(k),''',''asc_stress_Liu_V'');']);
	lon = lon-360;
	[lon,lat] = meshgrid(lon,lat);
	u = u';
	v = v';
	u(find(u==0))=nan;
	v(find(v==0))=nan;
	k
	fx = find(lon(1,:)>=-40 & lon(1,:)<=-32);
	fy = find(lat(:,1)>=-20 & lat(:,1)<=-10);
	
	lon = lon(fy,fx);
	lat = lat(fy,fx);
	u = u(fy,fx);
	v = v(fy,fx);
        % refazer programa de modo que o .mat tenha a data no nome do arquivo, sem ser em dias julianos !!!
        eval(['save tau0',num2str(k),'.mat lon lat u v'])
   
end

