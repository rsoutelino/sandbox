clear all;

rt=2;
t=0:1:2^(18)-1;
t=t*rt;
w1=2*pi/101/rt;
w2=2*pi/164/rt;
w3=2*pi/254/rt;
h1=.25*pi;
h2=0;
h3=pi;
A1=1;
A2=0.5;
A3=2;
ys1=A1*sin(w1*t+h1)+A2*sin(w2*t+h2)+A3*sin(w3*t+h3);
yn=30*randn(size(t));
yns1=yn+ys1;  
h1=h1+.5*pi;
h3=h3-.5*pi;
ys2=A3*sin(w1*t+h1)+A1*sin(w3*t+h3);
yn=30*randn(size(t));
yns2=yn+ys2;  

lt=length(t);

npf=35;

inxs=(npf-1)/2;

myfilter=blackman(npf);

myfilter=myfilter/sum(myfilter);

nrep=64;

nfft=round(lt/nrep);

fyb=zeros(nfft/2,nrep);

for j=1:nrep,

  i1=(j-1)*nfft+1;
  i2=i1+nfft-1;
  
  y1=conv(yns1(i1:i2),myfilter);
  y2=conv(yns2(i1:i2),myfilter);
  
  y1=y1(inxs+1:length(y1)-inxs);
  y2=y2(inxs+1:length(y2)-inxs);

  
  fy1=fft(y1); 
  fy2=fft(y2); 
  
  fy=2*fy2.*conj(fy1)/nfft.^2; 
  
  fyb(:,j)=fliplr(fy((length(fy)/2)+1:length(fy)))';

end

fyb=fyb';


ma=mean(sqrt(fyb.*conj(fyb)));
sa=std(sqrt(fyb.*conj(fyb)))/sqrt(nrep);

f=atan2(imag(fyb),real(fyb));
fu=cos(f);
fv=sin(f);
mf=atan2(mean(fv),mean(fu));


sf=atan2(std(imag(fyb)),std(real(fyb)))/sqrt(nrep);

freq=2*pi*[1:length(t(1:nfft))/2]/max(t(1:nfft));

mas=fliplr(sort(ma));
mas=mas(1:2);
idx=[find(ma==mas(1)) find(ma==mas(2))];
idm=idx-5;
idn=idx+5;

%portrait

subplot(2,1,1),loglog(freq,ma,'b',freq(idx),ma(idx),'or');
hold on
loglog([freq(idx); freq(idx)],[ma(idx)-sa(idx);ma(idx)+sa(idx)],'+r');
loglog(freq(idm), ma(idm),'om');
loglog([freq(idm); freq(idm)],[ma(idm)-sa(idm);ma(idm)+sa(idm)],'+m');
loglog(freq(idn), ma(idn),'om');
loglog([freq(idn); freq(idn)],[ma(idn)-sa(idn);ma(idn)+sa(idn)],'+m');
hold off
grid on;
axis([5e-3 5e-2 .2 2])
set(gca,'xtick',linspace(5e-3,5e-2,10))
set(gca,'ytick',linspace(.2,2,10))

ylabel('Potencia')
xlabel('Frequencia')
title(['Densidade de Potencia Espectral Cruzada (nrep=' num2str(nrep) ', npf=' ...
       num2str(npf) ', filtro blackman)' ])
text(0.0311,1.5e3,' <--0.0311')
text(0.0124,2.25e3,'0.0124--> ','HorizontalAlignment','right')




subplot(2,1,2),semilogx(freq,(180/pi)*mf,'b',freq(idx),(180/pi)*mf(idx),'or');
hold on
loglog([freq(idx); freq(idx)],(180/pi)*[mf(idx)-sf(idx);mf(idx)+sf(idx)],'+r');
hold off
grid on;
axis([5e-3 5e-2 -180 180])
text(0.0311,-90,' <-- -90','fontweight','bold')
text(0.0124,90,'90 --> ','HorizontalAlignment','right','fontweight','bold')
ylabel('Diferenca de fase')
xlabel('Frequencia')
set(gca,'xtick',linspace(5e-3,5e-2,10))


