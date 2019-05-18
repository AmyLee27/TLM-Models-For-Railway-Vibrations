function [M]=Melematrix1(x,ro)
% Calcula a matriz M do elemento
% - x - matriz [zi], com as coordenadas dos nós
% - ro - massa volúmica

% - Chamada de subrotinas
%   Call subroutines
Gauss1;
% - Alocação
%   Allocation
M=zeros(9,9);

% - Cálculo
%   Calculation
for ir=1:length(PG1(:,1));
    r=PG1(ir,2);
    [Jd]=Jacob1(x,r);
    [N]=Nmatrix1(r);
    Mi=N.'*N*ro*Jd;
    M=M+Mi;
end
