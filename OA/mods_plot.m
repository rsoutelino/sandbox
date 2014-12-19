close all;clear all;

lv=1; %level

[grdfile,grdpath] = uigetfile({'*.*'},' Pick a mod File ');
[header,hinfo,htype,z,temp,salt]=rhydro([grdpath,grdfile],1);

lat=hinfo(:,5);
lon=hinfo(:,4);


m_proj('mercator','lat',[-23 -5],'long',[-42 -30],'on');
m_gshhs_h('save','nwa_domain');

figure
hold on
m_plot(lon,lat,'*');
m_usercoast('nwa_domain','patch',[0.44 0.44 .44],'LineStyle','-');
m_grid('box','fancy','fontsize',10,'FontWeight','bold');
fn=[grdfile(1:16),'_map'];
%print( gcf, '-dpng', fn )
%print( gcf, '-depsc2', fn )

figure
hold on
m_usercoast('nwa_domain','patch',[0.44 0.44 .44],'LineStyle','-');
m_grid('box','fancy','fontsize',10,'FontWeight','bold');
[x,y]=m_ll2xy(lon,lat);
scatter(x,y,15,temp(:,lv),'filled')
colorbar

figure
hold on
m_usercoast('nwa_domain','patch',[0.44 0.44 .44],'LineStyle','-');
m_grid('box','fancy','fontsize',10,'FontWeight','bold');
[x,y]=m_ll2xy(lon,lat);
scatter(x,y,15,salt(:,lv),'filled')
colorbar
