%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      PROBLEMA TEORICO 2       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

H=4000;

n=0 % modo barotropico

Ai=sqrt((2)/(1-((sin(2*pi*n))/(2*pi*n))));

Fibt_0=Ai*cos(n*pi/4000*0)
Fibt_H=Ai*cos(n*pi/4000*(-4000))

n=1 % 1 modo baroclinico

Ai=sqrt((2)/(1-((sin(2*pi*n))/(2*pi*n))));
    
Fibc1_0=Ai*cos(n*pi/4000*0)
Fibc1_H=Ai*cos(n*pi/4000*(-4000))

n=2 % 2 modo baroclinico

Ai=sqrt((2)/(1-((sin(2*pi*n))/(2*pi*n))));
    
Fibc2_0=Ai*cos(n*pi/4000*0)
Fibc2_H=Ai*cos(n*pi/4000*(-4000))

n=3 % 3 modo baroclinico

Ai=sqrt((2)/(1-((sin(2*pi*n))/(2*pi*n))));
    
Fibc3_0=Ai*cos(n*pi/4000*0)
Fibc3_H=Ai*cos(n*pi/4000*(-4000))