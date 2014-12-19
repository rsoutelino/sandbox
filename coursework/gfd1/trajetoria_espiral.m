%%% Trajetoria da particula da lista 4 exerc. 5
x0=1;
y0=1;
r = -0.000001;
t = 1:10000000;
f0=-0.0001;

for i=1:length(t);
    x(i)=(x0/r)*(exp(r*t(i)))*(f0*sin(f0*t(i))+r*cos(f0*t(i)));
    y(i)=-(y0/f0)*(exp(r*t(i)))*(r*sin(f0*t(i))-f0*cos(f0*t(i)));
end
plot(x,y)