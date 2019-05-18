function [E G K lambda M Poisson Density Vs Vp Vr] = elastodynamic(I1,I2,I3,V1,V2,V3)
%VARIABLE NAMES - 'EE' 'G' 'K' 'l' 'M' 'v' 'p' 'Vs' 'Vp' 'Vr'

I = {I1;I2;I3};
V = [V1;V2;V3];

% EE - Young's modulus; 
% K - Elastic moduli (also called B: Bulk modulus);
% v - Possion's ratio v = 1/2 - E/6K; 
% p - density of the layer;
% G - Shear modulus; Vs = sqrt(G/p), which gives G = p*Vs^2
% l - lamda (lame's first parameter) lamda = 3*K*(3*K - EE)/(9*K - EE);
% M - P-wave (compressional)modulus Vp = sqrt(M/p),which gives M = p*Vp^2;
%     M = lamda + 2*G, based on the above eqns, M = 3*K*(3*K + EE)/(9*K - EE);
  
eqn = {    
'l = 3*K*(3*K-EE)/(9*K-EE)';...
'G = 3*K*EE/(9*K - EE)';...
'v = (3*K - EE)/(6*K)';...
'M = 3*K*(3*K +EE)/(9*K - EE)';...
'G = p*Vs^2';...
'M = p*Vp^2';...
'Vr = (0.87 + 1.12*v)*Vs/(1 + v)';
};
% Vr = (0.87 + 1.12*(l/(l+M)))*Vs/(1 + (l/(l+M)));
% or Vr = (0.87 + 1.12*v)*Vs/(1 + v)

for i = 1:3    
    eqnv = strcat(I{i},'=',num2str(V(i)));
    eqn = [eqn;eqnv];    
end

x= solve(eqn{1:numel(eqn)});
if isempty(x)
    disp('INCONSISTENT EQUATIONS');
    return;
end


E = (double(remmin(x.EE)));
G = (double(remmin(x.G)));
K = (double(remmin(x.K)));
lambda = (double(remmin(x.l)));
M = (double(remmin(x.M)));
Poisson = (double(remmin(x.v)));
Density = (double(remmin(x.p)));
Vs = (double(remmin(x.Vs)));
Vp = (double(remmin(x.Vp)));
Vr = (double(remmin(x.Vr)));

