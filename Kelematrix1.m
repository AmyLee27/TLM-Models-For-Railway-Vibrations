function [K]=Kelematrix1(x,k1,k2,E,niu,del,w)
% Calcula a matriz K do elemento
% - x - matriz [yi,zi], com as coordenadas dos nós
% - k - número de onda (rad/m)
% - E - módulo de Young;
% - niu - coeficiente de Poisson;
% - del - amortecimento;
% - w - frequencia (rad/s)

% - Chamada de subrotinas
[D]=Dmatrix(E,niu,del,w);
Gauss1;
% - Alocação
K=zeros(9,9);

% - Cálculo
for ir=1:length(PG1(:,1));
    r=PG1(ir,2);
    [B,Jd]=Bmatrix1(x,r,k1,k2);
    Ki=B'*D*B*Jd;
    K=K+Ki;
end
