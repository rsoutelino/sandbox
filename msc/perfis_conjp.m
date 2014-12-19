  clear;close all;clc;
  load mbar.mat
  alpha = mbar_model(1.40518e-4 ,lat,lon,bfrqs,profsl');
  clear lat lon 
  pp=[150 profsl'];
  
  load ocean_colormap
  load posicao.dat; 
  load monot2.mat
  lat=posicao(1:end-1,1) + posicao(1:end-1,2)/60;  
  lon=posicao(1:end-1,3) + posicao(1:end-1,4)/60; 
  profsl=posicao(1:end-1,5);
  [dst,ang]=sw_dist(lat,lon,'km');
  dist=[0 cumsum(dst)' 225];

  xx = 0:10:230;
  zz = -5000:100:0;
  [x,z] = meshgrid(xx,zz);
  contourf(x,z,z,[-5000:100:0]);colormap(ocean_colormap);hold on;shading flat;caxis([-6000 2000]);
  
  hold on
  fill([dist 0],[-150 -profsl' -profsl(end)],[.5 0.5 0.5]);
  
  hl1 = line(dist,-pp,'Color','k','linewidth',2);
  xlabel('Distancia [km]','fontsize',12)
  ylabel('Profundidade [m]','fontsize',12)
  axis([0 dist(end) pp(end)*-1 50])
  ax1 = gca;
  set(ax1,'XColor','k','YColor','k')
  

  ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','bottom','YAxisLocation','right','Color','none','XColor','k','YColor','r');
  hl2 = line(dist(2:end-1),alpha,'Color','r','linewidth',2,'Parent',ax2);
  set(gca,'XTickLabel',[]);
  ylabel('Potencial de Geracao','fontsize',12)
  axis([0 dist(end) -1 2.5])
  hl = line([0 dist(end)],[1 1],'Linestyle','--','Color','y','linewidth',1,'Parent',ax2);
  
  text([34.5],[-20],['Estacao 2'],'fontsize',12);

%  print(1,'-dpng',['figuras/perfil_tempsal_wb',XX]);

  pause
%  close all
  
end
