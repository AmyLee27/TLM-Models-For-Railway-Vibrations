function [K,M]=KMelematrix1(x,k1,k2,E,niu,del,ro,w)
% Calcula a matriz K do elemento
% - x - matriz [yi,zi], com as coordenadas dos n�s
% - k - n�mero de onda (rad/m)
% - E - m�dulo de Young;
% - niu - coeficiente de Poisson;
% - del - amortecimento;
% - w - frequencia (rad/s)
%ro - massa vol�mica;

% - Chamada de subrotinas
[D]=Dmatrix(E,niu,del,w);
Gauss1;
% - Aloca��o
K=zeros(9,9);
M=zeros(9,9);
% - C�lculo
for ir=1:length(PG1(:,1));
    r=PG1(ir,2);
    [B,Jd,N]=Bmatrix1(x,r,k1,k2);
    Ki=B'*D*B*Jd;
    K=K+Ki;
    Mi=N'*N*ro*Jd;
    M=M+Mi;
end
