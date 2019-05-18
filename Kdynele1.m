function [Kdy]=Kdynele1(x,prop,k1,k2,w)
% Calcula a matriz Kdyn do elemento
% - x - matriz [zi], com as coordenadas dos nós
% - ro - massa volúmica
% - k - número de onda (rad/m)
% - E - módulo de Young;
% - niu - coeficiente de Poisson;
% - del - amortecimento;
% - w - frequencia (rad/s)
% - prop- vector (linha, com as propriedades do elemento)
E=prop(1);
niu=prop(2);
del=prop(3);
ro=prop(4);

% - Chamada de subrotinas
[K,M]=KMelematrix1(x,k1,k2,E,niu,del,ro,w);
%[M]=Melematrix1(x,ro);

% - Cálculo
Kdy=K-w^2*M;
