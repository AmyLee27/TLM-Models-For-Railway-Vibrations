function [dNe]=dNelements1(r)
% Calcula os valores das funções dNi/dr, no ponto de gauss de coordenadas r
% Calculates the function values dNi/dr, in gauss point coordinates r
    
  dNe=[r-0.5;-2*r;r+0.5];