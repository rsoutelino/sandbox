clear all;close all;clc

eta0 = 1;
k = 0.0682;
H = 100;

for z = 0:0.5:H
    axis([-5 5 -6 1])
    axis('equal')
    eval(['ezplot(''x^2 + (z+',num2str(z),')^2 = exp(2*0.682*',num2str(-z),')'')']);hold on
%      eval(['ezplot(''x^2 + z^2 = exp(2*0.0682*',num2str(-z),')'')']);hold on
end


axis([-5 5 -6 1.2])
set(gca,'yticklabel',[])
set(gca,'xticklabel',[])
%  set(gca,'title',[])

pause

print -depsc exe4_fig1.eps
!epstopdf exe4_fig1.eps
!rm -rf exe4_fig1.eps

figure

for z = [0:0.8:3.9]
%      axis([-5 5 -6 1])
%      axis('equal')
    eval(['ezplot('' x^2/4 + (z+',num2str(z),')^2/(4 - ',num2str(z),') = 1 '')']);hold on
%      eval(['ezplot(''x^2 + z^2 = exp(2*0.0682*',num2str(-z),')'')']);hold on
end

axis([-3 3 -5 3])
set(gca,'yticklabel',[])
set(gca,'xticklabel',[])

pause

print -depsc exe4_fig2.eps
!epstopdf exe4_fig2.eps
!rm -rf exe4_fig2.eps

