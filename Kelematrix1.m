function [K]=Kelematrix1(x,k1,k2,E,niu,del,w)
% Calcula a matriz K do elemento
% - x - matriz [yi,zi], com as coordenadas dos n�s
% - k - n�mero de onda (rad/m)
% - E - m�dulo de Young;
% - niu - coeficiente de Poisson;
% - del - amortecimento;
% - w - frequencia (rad/s)

% - Chamada de subrotinas
[D]=Dmatrix(E,niu,del,w);
Gauss1;
% - Aloca��o
K=zeros(9,9);

% - C�lculo
for ir=1:length(PG1(:,1));
    r=PG1(ir,2);
    [B,Jd]=Bmatrix1(x,r,k1,k2);
    Ki=B'*D*B*Jd;
    K=K+Ki;
end
