
A=1;

dt=10;

t=-5*365:dt:5*365;

r=A*randn(size(t));

b=(20+t/1000)/1.058;

s=1.333*sin(2*pi*t/182.5)+0.5*sin(2*pi*t/365);

sst=r+b+s;

plot(t,sst,t,s+b)

vr=std(r)^2;

vb=std(b)^2;

vs=std(s)^2;

vsst=std(sst)^2;


[vr+vb+vs vsst]


sr=std(r);

sb=std(b);

ss=std(s);

ssst=std(sst);


[sr+sb+ss ssst]

xx=get(gca,'xlim');
yy=get(gca,'ylim');

text(xx(2)-1000,yy(1)+1,['Soma das variancias: ' num2str(vr+vb+vs)])
text(xx(2)-1000,yy(1)+2,['Variancia da soma: ' num2str(vsst)])
text(xx(1)+50,yy(2)-1,['Soma dos desvios: ' num2str(sr+sb+ss)])
text(xx(1)+50,yy(2)-2,['Desvio da soma: ' num2str(ssst)])

title('Desvio padrao X Variancia')
xlabel('Tempo(d)')
ylabel('Tempo(^oC)')
% print -depsc2 exe_pra_var.eps
