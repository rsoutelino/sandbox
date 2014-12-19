function err=vmode(rd);
global n2f2;
global po;
global f;
n=length(n2f2);
d=zeros(1,n);
f=ones(1,n);
nn=(n2f2(1:n-1)+n2f2(2:n))/2*.01;
dd=0;ff=1;g=1/rd/rd*.01;
for i=1:n-1
  dd=dd-g*ff;
  ff=ff+nn(i)*dd;
  d(i+1)=dd;
  f(i+1)=ff;
  end;
f=f(1:n-1);
f=f/sqrt(mean(f.*f));
xi0=mean(f.*f.*f);
%plot(f,-po);
err=dd;
end;

