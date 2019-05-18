function [K0,K1,K2,K3,K4,K5,M]=KMelematrixdiv1(x,prop)
% Calcula a matriz K do elemento
% - x - matriz [yi,zi], com as coordenadas dos nós
% - E - módulo de Young;
% - niu - coeficiente de Poisson;
% - del - amortecimento; % damping
% - w - frequencia (rad/s)
%ro - massa volúmica;

%%%%location of variables%%%
E=prop(1);
niu=prop(2);
del=prop(3);
ro=prop(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% - Chamada de subrotinas
[D]=Dmatrix(E,niu,del,1);
Gauss1;
% - Alocação
K0=zeros(9,9);
K1=zeros(9,9);
K2=zeros(9,9);
K3=zeros(9,9);
K4=zeros(9,9);
K5=zeros(9,9);
M=zeros(9,9);
% - Cálculo
for ir=1:length(PG1(:,1));
    r=PG1(ir,2);
    [B1,B2,B3,Jd,N]=Bmatrixdiv1(x,r);
    K0i=B1'*D*B1*Jd;
    K1i=B1'*D*B2*Jd-B2'*D*B1*Jd;
    K2i=B1'*D*B3*Jd-B3'*D*B1*Jd;
    K3i=B2'*D*B2*Jd;
    K4i=B3'*D*B3*Jd;
    K5i=B2'*D*B3*Jd+B3'*D*B2*Jd;
    K0=K0+K0i;
    K1=K1+K1i;
    K2=K2+K2i;
    K3=K3+K3i;
    K4=K4+K4i;
    K5=K5+K5i;
    
    Mi=N'*N*ro*Jd;
    M=M+Mi;
end
