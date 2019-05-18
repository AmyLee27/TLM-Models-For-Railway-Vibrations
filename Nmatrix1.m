function [N,Ne]=Nmatrix1(r)
% Calcula os valores das funções Ni, no ponto de gauss de cooerdenadas r
%Calculates the values of Ni functions in gauss point cooerdenadas r
  % r is the abscissar of Guass quadrature 
[Ne]=[0.5*r^2-0.5*r;1-r^2;0.5*r^2+0.5*r]; %shape functions expressed in LOCAL coordinates cf.in'introduction to FEM' Page 26
    
     [N]=[Ne(1),0,0, Ne(2),0,0,Ne(3),0,0;
          0,Ne(1),0,0,Ne(2),0,0,Ne(3),0,;
          0,0,Ne(1),0,0,Ne(2),0,0,Ne(3)];
      
          