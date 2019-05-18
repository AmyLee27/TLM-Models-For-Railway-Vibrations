function [Me,K0e,K1e,K2e,K3e,K4e,C0e]=beam1(E,nu,rho,A,Ixx,Izz,dofs)

mu=E/2/(1+nu);

Me=[rho*A  0  0;
    0  rho*A  0;
    0  0  rho*A];
K0e=zeros(3,3);
K1e=zeros(3,3);
K2e=[0  0  0;
     0 -E*A 0;
     0  0  0];
K3e=zeros(3,3);
K4e=[E*Izz 0  0;
     0   0   0;
     0   0  E*Ixx];
C0e=zeros(3,3);
switch dofs
    case 1
        DOF = [3];
    case 2
        DOF = [2,3];
    case 3
        DOF = [1,2,3];
end
Me=Me(DOF,DOF);
K0e=K0e(DOF,DOF);
K1e=K1e(DOF,DOF);
K2e=K2e(DOF,DOF);
K3e=K3e(DOF,DOF);
K4e=K4e(DOF,DOF);
C0e=C0e(DOF,DOF);