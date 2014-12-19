function [xi,yi,phii]=mirror_sc(xo,yo,phio);

% function mirror computes the image of a
% scalar quantity for OA  geostrophic streamline anlysis


% computing image position

m1=-(-45.5+53)/7.5;
b1= 7.5 -m1*(-53); 

m2=-1/m1;
b2= yo - m2*xo;

r= sqrt( m1*m1 + 1);

d= (-m1*xo + yo - b1)/r;

b3= -d*r +b1;

xi=-(b2-b3)/(m2-m1); 
yi= b3 + m1*xi;

%xi=xi-53;

% computing image scalar
phii=-phio;






