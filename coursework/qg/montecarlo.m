function [numlam,lamV,mlamc]=montecarlo(M,lamV,jmax);

%*************************************************
% Program to compute monte carlo simulations
% for EOF analysis
%
% it aims to determine which modes are
% statistically significant at the 95% level
%
% [numlam,lamV,mlamc]=montecarlo(M,jmax);
% M = Matrix of Data used in eof analysis
% jmax = number of simulations
% lamV = evalue of calculated eofs



% ICS Feb 2006
%************************************************

%%%%%%%%%%%%%%%% METHOD 1 %%%%%%%%%%%%%%%%%%%%%%

% method considering all data as statistically independent
% therefore, the size of random matriz is equal to the size
% of the data matriz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


V1=M;

[l1,l2]=size(V1);

lamc=zeros(l1,jmax);


for j=1:jmax,

	MC=randn(size(V1));  % use the randn function to assure a normal distr.

[lam,F,A]=eoft(MC); % ***COMPUTE EOFs HERE***

lamc(:,j)=lam;

end

mlamc=mean(lamc,2);
dlamc=std(lamc,0,2);

numlam=1:10;

figure
     plot(numlam, lamV(end:-1:end-9)*100,'b-*')
hold on
     plot(numlam, mlamc(end:-1:end-9)*100,'rd--')
     plot(numlam, lamV(end:-1:end-9)*100+2*dlamc(end:-1:end-9)*100,'c:')
     plot(numlam, lamV(end:-1:end-9)*100-2*dlamc(end:-1:end-9)*100,'c:')
hold off
xlabel('Numero do Modo EOF')
ylabel('Percentual da Variancia Explicada')
legend('Real','Simulado')
     title('Metodo 1 - Serie Simulada do Tamanho da Serie Real') 
     axis([1 10 0 100])
%pause


%%%%%%%%%%%%%%%% METHOD 2 %%%%%%%%%%%%%%%%%%%%%%

% method considering only the statistically independent data
% in the series.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1=effss(V1); % calculate no. of statistically independent values
n2=ceil(n1/10);


for j=1:jmax,

	MC2=randn(10,n2);  % use the randn function to assure a normal distr.

[lam,F,A]=eoft(MC); % ***COMPUTE EOFs HERE***

lamc2(:,j)=lam;

end

mlamc2=mean(lamc2,2);
dlamc2=std(lamc2,0,2);

figure
     plot(numlam, lamV(end:-1:end-9)*100,'b-*')
hold on
     plot(numlam, mlamc2(end:-1:end-9)*100,'rd--')
     plot(numlam, lamV(end:-1:end-9)*100+2*dlamc2(end:-1:end-9)*100,'c:')
     plot(numlam, lamV(end:-1:end-9)*100-2*dlamc2(end:-1:end-9)*100,'c:')
hold off
xlabel('Numero do Modo EOF')
ylabel('Percentual da Variancia Explicada')
legend('Real','Simulado')
title('Metodo 1 - Serie Simulada com Graus de Liberdade Calculados da Serie Real') 
     axis([1 10 0 100])
