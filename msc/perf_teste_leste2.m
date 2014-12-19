% plota perfis hidrograficos para visualizacao rapida
% Oceano Leste 2 - Mestrado

clear all;close all;clc

for i=[1:4 7:102 105:112]
      
      eval(['load ../../../dados/leste2/ctd/filtrados/lesteII_ctd',num2str(i),'.dat;'])
      eval(['ctd = lesteII_ctd',num2str(i),';'])
      eval(['clear lesteII_ctd',num2str(i)])

      p = ctd(:,1);
      t = ctd(:,2);
      s = ctd(:,3);
      
figure(1)
   subplot(121)
      plot(t,-p)
      tit=['Estacao',num2str(i)];
      title(tit)
   subplot(122)
      plot(s,-p,'r')
pause
end
