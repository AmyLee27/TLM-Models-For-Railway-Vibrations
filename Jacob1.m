function [Jd,dNx]=Jacob1(x,r)
% Calcula o determinante do Jacobiano e a derivada em ordem as coordenadas
% naturais das funções de forma

% Calculates the determinant of the Jacobian and the derivative in order to
% coordinate natural form of functions

[dNe]=dNelements1(r);
J=x.'*dNe; % the Jacobian matrix
Jd=det(J); % the determinant of the jacobian matrix
dNx=dNe*inv(J); % cf. in 'introduction to FEM' Page 24 eqn(4.17)
% The derivatives in respect to the GLOBAL coordinates 