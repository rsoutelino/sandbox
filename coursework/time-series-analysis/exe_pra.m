clear all;

if 1==1,
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
save exe_pra_1 yns t rt
end

clear all;

load exe_pra_1

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
  fy=fy.*conj(fy); fy2=fy;
  fyb(:,j)=fliplr(fy((length(fy)/2)+1:length(fy)))';
end

fyb=fyb';

freq=2*pi*[1:length(t(1:nfft))/2]/max(t(1:nfft));

mf=mean(fyb);
sf=std(fyb);

subplot(2,1,1),loglog(freq,mf,'b');grid on;axis('tight')
subplot(2,1,2),loglog(freq,mf,'b',freq,mf+1*sf,'y',freq,mf+2*sf,'m',freq,mf+3*sf,'r');grid on;axis('tight')


