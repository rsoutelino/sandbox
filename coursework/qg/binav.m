function cfit=binav(pi,p,c);

n=length(pi)-1;
cfit=zeros(n,1);
for i=1:n
  cs=c(p>=pi(i) & p<pi(i+1) & ~isnan(c));
  cfit(i)=sum(cs)/sum(ones(size(cs)));
  end
return
