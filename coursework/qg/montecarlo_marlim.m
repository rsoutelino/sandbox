%*************************************************
% program to compute monte carlo simulations
% for the Marlim EOF analysis

% it aims to determine which modes are
% statistically significant at the 95% level

% ICS Feb 2006
%************************************************

%%%%%%%%%%%%%%%% METHOD 1 %%%%%%%%%%%%%%%%%%%%%%

% method considering all data as statistically independent
% therefore, the size of random matriz is equal to the size
% of the data matriz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load matriz_V

jmax=1000;

[l1,l2]=size(V1);

lamc=zeros(l1,jmax);


for j=1:jmax,

	MC=randn(size(V1));  % use the randn function to assure a normal distr.

[lam,F,A]=eoft(MC); % ***COMPUTE EOFs HERE***

lamc(:,j)=lam;
k=j-1;
disp(k)

end

mlamc=mean(lamc,2);
dlamc=std(lamc,0,2);

numlam=1:10;

figure(1)
     plot(numlam, lamV(10:-1:1)*100,'b*',numlam,lamV(10:-1:1)*100,'b-')
hold on
     plot(numlam, mlamc(10:-1:1)*100,'rd',numlam,mlamc(10:-1:1)*100,'r--')
     plot(numlam, lamV(10:-1:1)*100+2*dlamc(10:-1:1)*100,'c:')
     plot(numlam, lamV(10:-1:1)*100-2*dlamc(10:-1:1)*100,'c:')
hold off
xlabel('EOF Mode Number')
ylabel('Percentual Variance Explained')
     title('Method 1') 
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

figure(2)
     plot(numlam, lamV(10:-1:1)*100,'b*',numlam,lamV(10:-1:1)*100,'b-')
hold on
     plot(numlam, mlamc2(10:-1:1)*100,'rd',numlam,mlamc2(10:-1:1)*100,'r--')
     plot(numlam, lamV(10:-1:1)*100+2*dlamc2(10:-1:1)*100,'c:')
     plot(numlam, lamV(10:-1:1)*100-2*dlamc2(10:-1:1)*100,'c:')
hold off
xlabel('EOF Mode Number')
ylabel('Percentual Variance Explained') 
     title('Method 2')
     axis([1 10 0 100])


