function [B1,B2,B3,Jd,N]=Bmatrixdiv1(x,r)
% Calcula a matriz B do elemento
% - r - coordenada local
% - s coordenada local
% - x - matriz [yi,zi], com as coordenadas dos nós
% - k - número de onda (rad/m)
% Calculate the strain-displacement matrix, in which the derivatives are
% expressed in GLOBAL coordinates

[Jd,dNx]=Jacob1(x,r);
[N,Ne]=Nmatrix1(r);

%%% Matrix B0%%%%%%
B1=[0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0, 0,0,0;
    0,0,dNx(1),0,0,dNx(2),0,0,dNx(3);
    0,0,0,0,0,0,0,0,0;
    0,dNx(1),0, 0,dNx(2),0, 0,dNx(3),0;
    dNx(1),0,0,dNx(2),0,0,dNx(3),0,0];

%%%Matrix B1%%%%%
B2=[Ne(1),0,0,Ne(2),0,0,Ne(3),0,0;
    0, 0,0,0,0,0, 0,0,0;
    0,0,0,0,0,0,0,0,0;
    0,Ne(1),0,0,Ne(2),0,0,Ne(3),0;
    0,0,0, 0,0,0, 0,0,0;
    0,0,Ne(1),0,0,Ne(2),0,0,Ne(3)];

%%%%%Matrix B2%%%%%
B3=[0,0,0,0,0,0,0,0,0;
    0,Ne(1),0,0,Ne(2),0, 0,Ne(3),0;
    0,0,0,0,0,0,0,0,0;
    Ne(1),0,0,Ne(2),0,0,Ne(3),0,0;
    0,0,Ne(1),0,0,Ne(2), 0,0,Ne(3);
    0,0,0,0,0,0,0,0,0];
    