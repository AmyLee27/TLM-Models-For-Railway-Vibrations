function [D]=Dmatrix(E,niu,del,w)
% Calcula a matriz D

% calculate lamda and mu (Lames constants) and also elasticity matrix E(or D in this code)
% cf. in 'Introduction to the finite element method' Chapter 3 Page 16
% where E is elasticity modulus
       la=niu*E*(1+2*i*del*sign(w))/((1+niu)*(1-2*niu)); 
       mi=E*(1+2*i*del*sign(w))/(2*(1+niu));
       D=[la+2*mi,la,la,0,0, 0;la,la+2*mi,la,0,0,0;la,la,la+2*mi,0,0,0;0,0,0,mi,0,0;0,0,0,0,mi,0;0,0,0,0,0,mi]; 
          