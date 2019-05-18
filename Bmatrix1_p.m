function [B,Jd,N] = Bmatrix1_p(x,r,k1,k2) 
  
% Calcula a matriz B do elemento
% Calculates the B matrix element

% - r - coordenada local
%       coordinate location

% - s coordenada local
% - x - matriz [yi,zi], com as coordenadas dos nós
%                       with the coordinates of the nodes
                            
% - k - número de onda (rad/m) 
%       wave number                                                 

[Jd,dNx]=Jacob1(x,r);
[N,Ne]=Nmatrix1(r);


 B=[-1i*k1*Ne(1),0,0,-1i*k1*Ne(2),0,0,-1i*k1*Ne(3),0,0;
     0, -1i*k2*Ne(1),0,0, -1i*k2*Ne(2),0, 0, -1i*k2*Ne(3),0;
     0,0,dNx(1),0,0,dNx(2),0,0,dNx(3);
     -1i*k2*Ne(1), -1i*k1*Ne(1),0,-1i*k2*Ne(2), -1i*k1*Ne(2),0,-1i*k2*Ne(3), -1i*k1*Ne(3),0;
     0,dNx(1),-1i*k2*Ne(1), 0,dNx(2),-1i*k2*Ne(2), 0,dNx(3),-1i*k2*Ne(3);
     dNx(1),0,-1i*k1*Ne(1),dNx(2),0,-1i*k1*Ne(2),dNx(3),0,-1i*k1*Ne(3)];
    