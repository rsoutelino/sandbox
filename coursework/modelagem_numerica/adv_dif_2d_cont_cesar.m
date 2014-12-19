% Programa para adveccao e difusao 2D
% de vazamento continuo de poluente no mar

%%%%constantes do modelo
mmax=33;  % numero de colunas da grade
lmax=11;  % numero de linhas na grade
nmax=90;  % numero de instantes de tempo
dx=1;     % espacamento, em Km
dx2=dx*dx;
dt=1;     % passo de tempo, em horas
D=3.6e-5; % coeficiente de difusao, em km^2/horas
lvaz=6; mvaz=8; xvaz=180; % posicao e valor do vazamento

%%%%Campo de velocidades, em Km /h
u=ones(lmax,mmax)*.54; 
u(1,:)=0;
u(lmax,:)=0;
v=zeros(lmax,mmax);

%%%%%Condicoes iniciais
fadv=zeros(lmax,mmax);
fant=zeros(lmax,mmax);
fatu=zeros(lmax,mmax);
fren=zeros(lmax,mmax);
fvaz=zeros(lmax,mmax);
fvaz(lvaz,mvaz)=xvaz;
fant=fvaz;
fatu=fvaz;

%%%%%Equacao da adveccao - difusao 2D
for t=2:nmax
    fadv(2:end-1,2:end-1)=-((v(2:end-1,2:end-1).*(fatu(3:end,2:end-1)-fatu(1:end-2,2:end-1)))/(2*dx)...
        +u(2:end-1,2:end-1).*(fatu(2:end-1,3:end)-fatu(2:end-1,1:end-2))/(2*dx));
    as=find(fadv<0);
    fadv(as)=0;
    fren(2:end-1,2:end-1)=fant(2:end-1,2:end-1)+2*dt*...
    (D*((fant(3:end,2:end-1)+fant(1:end-2,2:end-1)-2*fant(2:end-1,2:end-1))/dx2+ ...
    (fant(2:end-1,3:end)+fant(2:end-1,1:end-2)-2*fant(2:end-1,2:end-1))/dx2)+...
    fadv(2:end-1,2:end-1));
    fren=fren+fvaz;
    %fren=fren+fvaz;
    contourf(fren)
    title(['Adveccao e Difusao 2D de vazamento continuo - concentracoes em t = ',...
          num2str(t),' horas'])
    xlabel('Km')
    ylabel('Km')
    colorbar
    pause
    fant=fatu;
    fatu=fren;
     
end

