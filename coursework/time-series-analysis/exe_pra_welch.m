clear all;

rt=2;
t=0:1:2^(18)-1;
t=t*rt;
w1=2*pi/101/rt;
w2=2*pi/104/rt;
w3=2*pi/254/rt;
h1=.25*pi;
h2=0;
h3=pi;
A1=1;
A2=1.5;
A3=2;
ys=A1*sin(w1*t+h1)+A2*sin(w2*t+h2)+A3*sin(w3*t+h3);
yn=20*randn(size(t));
yns=yn+ys;  

lt=length(t);

npf=51;

inxs=(npf-1)/2;

myfilter=blackman(npf);

myfilter=myfilter/sum(myfilter);

nrep=32;

nfft=round(lt/nrep);

fyb=zeros(nfft/2,nrep);

for j=1:nrep,

  i1=(j-1)*nfft+1;
  i2=i1+nfft-1;
  
  y=conv(yns(i1:i2),myfilter);
  
  y=y(inxs+1:length(y)-inxs);
  
  fy=fft(y); 
  Am=2*fy.*conj(fy)/(nfft.^2);
  Amb(:,j)=fliplr(Am((length(Am)/2)+1:length(Am)))';
end

Amb=Amb';

freq=2*pi*[1:length(t(1:nfft))/2]/max(t(1:nfft));

mf=mean(Amb);
sf=std(Amb);
%----------------------------------------------------
%portrait
aa=[ 1e2 2e4 NaN   1e2 2e4 NaN  1e2 2e4];
aw=[0.0311 0.0311 NaN 0.0302 0.0302 NaN 0.0124 0.0124];

subplot(2,1,1),loglog(freq,mf,'b');
hold on;plot(aw,aa,'--k');hold off
grid on;
axis([1e-2 4e-2 5e-2 3])
set(gca,'xtick',linspace(1e-2,4e-2,10))
ylabel('Potencia')
xlabel('Frequencia')
title(['Densidade de Potencia Espectral (nrep=' num2str(nrep) ', npf=' ...
       num2str(npf) ', filtro blackman)' ])
text(0.0311,1.1e4,'<--0.0311')
text(0.0302,1.1e4,'0.0302-->','HorizontalAlignment','right')
text(0.0124,1.1e4,'<--0.0124')



subplot(2,1,2),loglog(freq,mf,'b',freq,mf+1*sf,'y',freq,mf+2*sf,'m',freq, ...
                      mf+3*sf,'r');
hold on;plot(aw,aa,'--k');hold off
grid on;
axis([1e-2 5e-2 1e-2 3])
set(gca,'xtick',linspace(1e-2,4e-2,10))
ylabel('Potencia')
xlabel('Frequencia')
legend('Espectro','68%','95%','99.7%')
