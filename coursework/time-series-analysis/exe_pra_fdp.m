clear all
m=6;
s=4;

% 1000 repeticoes de 1120 amostras cada

Z=rand(1000,1120);
mZ=mean(Z);
subplot(2,2,1),hist(Z(123,:),100)
subplot(2,2,2),hist(Z(234,:),100)
subplot(2,2,3),hist(Z(666,:),100)
subplot(2,2,4),hist(mZ,100)

