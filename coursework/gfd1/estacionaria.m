%%%%%%% Exercicio 5 onda estacionaria %%%%%%%

N0=0.5;
w=0.43;
k=0.014;

for x=1:1:1000;
    Ne0(x)=2*N0*cos(k*x)*cos(-w*0);
    Ne2(x)=2*N0*cos(k*x)*cos(-w*2);
    Ne8(x)=2*N0*cos(k*x)*cos(-w*8);
    Ne9(x)=2*N0*cos(k*x)*cos(-w*9);
end

plot(Ne0,'k:')
hold
plot(Ne2,'b')
plot(Ne8,'k')
plot(Ne9,'b:')
title('Elevacao de uma onda estacionaria')
xlabel('Eixo x (m)')
ylabel('Eta (m)')
hold off
print -dtiff onda_est.tif