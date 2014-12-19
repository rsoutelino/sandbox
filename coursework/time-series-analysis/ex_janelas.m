%%%  plota exemplos de janelas moveis %%%%

myfilter1=triang(31);
myfilter2=hanning(31);
x=0:30;

figure(1)
set(gcf,...
        'Color',[1 1 1],...
        'InvertHardcopy','on',...
        'PaperUnits','inches',...
        'Units','inches',...
        'PaperOrientation','portrait',...
        'PaperPosition',[0 0 8.5 11],...
        'PaperPositionMode','manual',...
        'PaperType','usletter',...
        'Position',[.2 .2 8.5 11],...
        'ShareColors','off',...
        'Clipping','on');

subplot(2,1,1)
plot(x,myfilter1,'r','linewidth',2);
hold on
plot(15*ones(11),0:.1:1,'k','linewidth',0.5)
axis([-1 32 0 1.2])
title('Exemplo de janela triangular de 31 pontos','fontweight','bold','fontsize',12)
xlabel('Numero de Amostras')
ylabel('Pesos Atribuídos')

subplot(2,1,2)
plot(x,myfilter2,'linewidth',2);
hold on
plot(15*ones(11),0:.1:1,'k','linewidth',0.5)
title('Exemplo de janela gaussiana de 31 pontos','fontweight','bold','fontsize',12)
xlabel('Numero de Amostras')
ylabel('Pesos Atribuídos')
axis([-1 32 0 1.2])

print -depsc ex_janelas.eps
!epstopdf ex_janelas.eps