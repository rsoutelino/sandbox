%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ROTINA DE TESTE DA ANALISE OBJETIVA, USANDO SCALOA %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

%%% carregando os dados: distancia, profundidade, temperatura (secao vertical)
load exe.mat;
T=T';
dist=dist';
prof=prof';

%%% criando grade para qual os dados serao interpolados

x=0:.5:max(dist);
z=max(prof):-1:min(prof);

lx=length(x);
lz=length(z);

[xx,zz]=meshgrid(x,z);
xx=reshape(xx,1,lx*lz);
zz=reshape(zz,1,lx*lz);

corrlen=5%input('Entre com o comprimento de correlacao')
err=0.07%input('Entre com a variancia do erro aleatorio')

[To,Eo]=scaloa(xx,zz,dist,prof,T,corrlen,err);

%%% reshapeando pra plotar

To=reshape(To,lz,lx);
Eo=reshape(Eo,lz,lx);
xx=reshape(xx,lz,lx);
zz=reshape(zz,lz,lx);

dist=reshape(dist,14,4);
prof=reshape(prof,14,4);
T=reshape(T,14,4);

%%% plotando e comparando

lT=22:0.05:28.5;
lE=0.001:0.0005:0.0517;

subplot(2,2,1)
contourf(dist,prof,T,lT);shading flat;hold on;colorbar
xlabel('Distancia')
ylabel('Profundidade')
title('Dados brutos')

subplot(2,2,2)
contourf(xx,zz,To,lT);shading flat;hold on;colorbar
xlabel('Distancia')
ylabel('Profundidade')
title('Interpolacao Objetiva')

subplot(2,2,4)
contourf(xx,zz,Eo,lE);shading flat;hold on;colorbar
xlabel('Distancia')
ylabel('Profundidade')
title('Erro')


