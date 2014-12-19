clear all;
mt=10000;
t=0:1:mt;
w1=2*pi/(mt/20);
w2=2*pi/(mt/25);
w3=2*pi/(mt*2);
y1=sin(w1*t);
y2=sin(w2*t);
y3=sin(w3*t);
no=(sqrt(2)/4)*randn(size(t));

y=2+y3.*(y1+y2+no);

m=5;g=gray(m);g=g(1:m-1,:);g=flipud(g)

figure(1)
set(gcf,'DefaultAxesColororder',g);
ta=t(1980:1990);
ya=y(1980:1990);
p=polyfit(ta,ya,1);
yla=polyval(p,ta);
plot(ta,ya,ta,yla,'--');axis('tight')
xlabel('Tempo');ylabel('Vento')
title('Tendencia linear?')
print -depsc2 exe_pra_len1.eps

figure(2)
set(gcf,'DefaultAxesColororder',g);
tb=t(1900:2100);
yb=y(1900:2100);
p=polyfit(tb,yb,1);
ylb=polyval(p,tb);
plot(tb,yb,tb,ylb,'--',ta,ya,ta,yla,'--');
axis('tight');
xlabel('Tempo');ylabel('Vento')
title('Tendencia linear.')
print -depsc2 exe_pra_len2.eps

figure(3)
set(gcf,'DefaultAxesColororder',g);
tc=t(1700:2500);
yc=y(1700:2500);
ylc=2+1.2*sin(-1.5+tc*2*pi/470); % ajustado por tentativa e erro
plot(tc,yc,tc,ylc,'--',tb,yb,tb,ylb,'--')
axis('tight');
xlabel('Tempo');ylabel('Vento')
title('Era uma onda...')
print -depsc2 exe_pra_len3.eps

figure(4)
set(gcf,'DefaultAxesColororder',g);
td=t(1000:3000);
yd=y(1000:3000);
plot(td,yd,tb,yb,tb,ylb,'--')
axis('tight');
xlabel('Tempo');ylabel('Vento')
title('Duas ondas? Batimento?')
print -depsc2 exe_pra_len4.eps

figure(5)
set(gcf,'DefaultAxesColororder',g);
plot(t,y,td,yd,tc,yc,tb,yb,ta,ya)
axis('tight');
xlabel('Tempo');ylabel('Vento')
title('Serie curta, conclusoes erradas')
print -depsc2 exe_pra_len5.eps

