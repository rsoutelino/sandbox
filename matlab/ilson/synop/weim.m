function [y,w]=weim(N,wt,x)

% function [y,w]=weim(N,wt,x)
%
% it computes a weighted average using the window function
% this filter is designed for odd weight numbers only
% N is numbers of weigths

w=window(N,wt);

r=size(x);
if r(1) > 1, x=x'; end

ln=(N-1)/2;
lx=length(x);
lf=lx-ln+1;
y=zeros(size(x));

for i=1:lx,

  if i <= ln,
  
   y(i)=sum(x(1:ln+i).*w(ln+2-i:N))/sum(w(ln+2-i:N));
  
  elseif ((i > ln) & (i < lf)),
  
   y(i)=sum(x(i-ln:i+ln).*w)/sum(w);

  else  % i >=lf 

   y(i)=sum(x(i-ln:lx).*w(1:length(i-ln:lx)))/sum(w(1:length(i-ln:lx)));

  
  end
end
