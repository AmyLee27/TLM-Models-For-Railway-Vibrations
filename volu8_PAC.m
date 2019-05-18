function [Me,K0e,K1e,K2e,K3e,K4e,C0e]=volu8_PAC(x,E,nu,rho,dofs)

lambda=nu*E/(1+nu)/(1-2*nu);
mu=E/2/(1+nu);

C=[lambda+2*mu lambda      lambda      0  0  0
    lambda      lambda+2*mu lambda      0  0  0
    lambda      lambda      lambda+2*mu 0  0  0
    0           0           0           mu 0  0
    0           0           0           0  mu 0
    0           0           0           0  0  mu];

nXi1D=10;
nXi=nXi1D^2;
xi1D=[0.97390652851717;0.86506336668898;0.67940956829902;0.43339539412925;0.14887433898163;-0.14887433898163;-0.43339539412925;-0.67940956829902;-0.86506336668898;-0.97390652851717];
H1D =[0.06667134430869;0.14945134915058;0.21908636251598;0.26926671931000;0.29552422471475; 0.29552422471475; 0.26926671931000; 0.21908636251598; 0.14945134915058; 0.06667134430869];
xi=zeros(nXi,2);
H=zeros(nXi,1);
for iXi=1:nXi1D
    for jXi=1:nXi1D
        xi(nXi1D*(iXi-1)+jXi,1)= xi1D(iXi);
        xi(nXi1D*(iXi-1)+jXi,2)= xi1D(jXi);
        H(nXi1D*(iXi-1)+jXi)=H1D(iXi)*H1D(jXi);
    end
end

N=zeros(nXi,8);
% N(:,1)=0.25*(1.0-xi(:,1)).*(xi(:,2)-1.0).*(xi(:,1)+xi(:,2)+1.0); %1
% N(:,2)=0.5*(1.0-xi(:,2)).*(1.0-xi(:,1).^2); %2
% N(:,3)=0.25*(1.0+xi(:,1)).*(xi(:,2)-1.0).*(-xi(:,1)+xi(:,2)+1.0); %3
% N(:,4)=0.5*(1.0+xi(:,1)).*(1.0-xi(:,2).^2); %4
% N(:,5)=0.25*(1.0+xi(:,1)).*(1.0+xi(:,2)).*(xi(:,1)+xi(:,2)-1.0); %5
% N(:,6)=0.5*(1.0+xi(:,2)).*(1.0-xi(:,1).^2); %6
% N(:,7)=-0.25*(1.0-xi(:,1)).*(1.0+xi(:,2)).*(xi(:,1)-xi(:,2)+1.0); %7
% N(:,8)=0.5*(1.0-xi(:,1)).*(1.0-xi(:,2).^2); %8

N(:,1)=0.25*(1.0-xi(:,1)).*(xi(:,2)-1.0).*(xi(:,1)+xi(:,2)+1.0); %1
N(:,5)=0.5*(1.0-xi(:,2)).*(1.0-xi(:,1).^2); %5
N(:,2)=0.25*(1.0+xi(:,1)).*(xi(:,2)-1.0).*(-xi(:,1)+xi(:,2)+1.0); %2
N(:,6)=0.5*(1.0+xi(:,1)).*(1.0-xi(:,2).^2); %6
N(:,3)=0.25*(1.0+xi(:,1)).*(1.0+xi(:,2)).*(xi(:,1)+xi(:,2)-1.0); %3
N(:,7)=0.5*(1.0+xi(:,2)).*(1.0-xi(:,1).^2); %7
N(:,4)=-0.25*(1.0-xi(:,1)).*(1.0+xi(:,2)).*(xi(:,1)-xi(:,2)+1.0); %4
N(:,8)=0.5*(1.0-xi(:,1)).*(1.0-xi(:,2).^2);  %8


dN=zeros(nXi,16);
% dN(:,1)=0.25*(1.0-xi(:,2)).*(2*xi(:,1)+xi(:,2));
% dN(:,9)=0.25*(1.0-xi(:,1)).*(2*xi(:,2)+xi(:,1));
% dN(:,2)=-1.0*xi(:,1).*(1.0-xi(:,2));
% dN(:,10)=-0.5*(1.0-xi(:,1).^2);
% dN(:,3)= -0.25*(1.0-xi(:,2)).*(xi(:,2)-2*xi(:,1));
% dN(:,11)=0.25*(1.0+xi(:,1)).*(2*xi(:,2)-xi(:,1));
% dN(:,4)=0.5*(1.0-xi(:,2).^2);
% dN(:,12)=-1.0*xi(:,2).*(1.0+xi(:,1));
% dN(:,5)= 0.25*(1.0+xi(:,2)).*(2*xi(:,1)+xi(:,2));
% dN(:,13)=0.25*(1.0+xi(:,1)).*(2*xi(:,2)+xi(:,1));
% dN(:,6)=-1.0*xi(:,1).*(1.0+xi(:,2));
% dN(:,14)= 0.5*(1.0-xi(:,1).^2);
% dN(:,7)=-0.25*(1.0+xi(:,2)).*(xi(:,2)-2*xi(:,1));
% dN(:,15)=0.25*(1.0-xi(:,1)).*(2*xi(:,2)-xi(:,1));
% dN(:,8)=-0.5*(1.0-xi(:,2).^2);
% dN(:,16)=-1.0*xi(:,2).*(1.0-xi(:,1));

dN(:,1)=0.25*(1.0-xi(:,2)).*(2*xi(:,1)+xi(:,2));
dN(:,9)=0.25*(1.0-xi(:,1)).*(2*xi(:,2)+xi(:,1));
dN(:,5)=-1.0*xi(:,1).*(1.0-xi(:,2));
dN(:,13)=-0.5*(1.0-xi(:,1).^2);
dN(:,2)= -0.25*(1.0-xi(:,2)).*(xi(:,2)-2*xi(:,1));
dN(:,10)=0.25*(1.0+xi(:,1)).*(2*xi(:,2)-xi(:,1));
dN(:,6)=0.5*(1.0-xi(:,2).^2);
dN(:,14)=-1.0*xi(:,2).*(1.0+xi(:,1));
dN(:,3)= 0.25*(1.0+xi(:,2)).*(2*xi(:,1)+xi(:,2));
dN(:,11)=0.25*(1.0+xi(:,1)).*(2*xi(:,2)+xi(:,1));
dN(:,7)=-1.0*xi(:,1).*(1.0+xi(:,2));
dN(:,15)= 0.5*(1.0-xi(:,1).^2);
dN(:,4)=-0.25*(1.0+xi(:,2)).*(xi(:,2)-2*xi(:,1));
dN(:,12)=0.25*(1.0-xi(:,1)).*(2*xi(:,2)-xi(:,1));
dN(:,8)=-0.5*(1.0-xi(:,2).^2);
dN(:,16)=-1.0*xi(:,2).*(1.0-xi(:,1));


Jac=zeros(nXi,2,2);
for iXi=1:nXi
    for iNod=1:8
        Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)*x(iNod,1);
        Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)*x(iNod,2);
        Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,8+iNod)*x(iNod,1);
        Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,8+iNod)*x(iNod,2);
    end
end

Me =zeros(24,24);
K0e=zeros(24,24);
K1e=zeros(24,24);
K2e=zeros(24,24);
K3e=zeros(24,24);
K4e=zeros(24,24);
C0e=zeros(24,24);
for iXi=1:nXi
    JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2);Jac(iXi,2,1),Jac(iXi,2,2)];
    detJac=det(JacUtil);
    DN1= JacUtil\[dN(iXi,1);dN(iXi,9)];
    DN2= JacUtil\[dN(iXi,2);dN(iXi,10)];
    DN3= JacUtil\[dN(iXi,3);dN(iXi,11)];
    DN4= JacUtil\[dN(iXi,4);dN(iXi,12)];
    DN5= JacUtil\[dN(iXi,5);dN(iXi,13)];
    DN6= JacUtil\[dN(iXi,6);dN(iXi,14)];
    DN7= JacUtil\[dN(iXi,7);dN(iXi,15)];
    DN8= JacUtil\[dN(iXi,8);dN(iXi,16)];
      
    B1=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,DN1(1),0,0,DN2(1),0,0,DN3(1),0,0,DN4(1),0,0,DN5(1),0,0,DN6(1),0,0,DN7(1),0,0,DN8(1),0;
    0,0,DN1(2),0,0,DN2(2),0,0,DN3(2),0,0,DN4(2),0,0,DN5(2),0,0,DN6(2),0,0,DN7(2),0,0,DN8(2);
    DN1(1),0,0,DN2(1),0,0,DN3(1),0,0,DN4(1),0,0,DN5(1),0,0,DN6(1),0,0,DN7(1),0,0,DN8(1),0,0;
    0,DN1(2),DN1(1), 0,DN2(2),DN2(1),0,DN3(2),DN3(1),0,DN4(2),DN4(1),0,DN5(2),DN5(1),0,DN6(2),DN6(1),0,DN7(2),DN7(1),0,DN8(2),DN8(1);
    DN1(2),0,0,DN2(2),0,0,DN3(2),0,0,DN4(2),0,0,DN5(2),0,0,DN6(2),0,0,DN7(2),0,0,DN8(2),0,0];

    B2=[N(iXi,1),0,0,N(iXi,2),0,0,N(iXi,3),0,0,N(iXi,4),0,0,N(iXi,5),0,0,N(iXi,6),0,0,N(iXi,7),0,0,N(iXi,8),0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,N(iXi,1),0,0,N(iXi,2),0,0,N(iXi,3),0,0,N(iXi,4),0,0,N(iXi,5),0,0,N(iXi,6),0,0,N(iXi,7),0,0,N(iXi,8),0;
    0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,N(iXi,1),0,0,N(iXi,2),0,0,N(iXi,3),0,0,N(iXi,4),0,0,N(iXi,5),0,0,N(iXi,6),0,0,N(iXi,7),0,0,N(iXi,8)];
    
    Nxi_1=[N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0        0
        0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4) 0
        0        0        N(iXi,1) 0        0        N(iXi,2) 0        0        N(iXi,3) 0        0        N(iXi,4)];
    
    Nxi_2=[N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0        0
        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8) 0
        0        0        N(iXi,5) 0        0        N(iXi,6) 0        0        N(iXi,7) 0        0        N(iXi,8)];
    
    Nxi = [Nxi_1 Nxi_2];
    
    Me= Me +H(iXi)*detJac*(rho*(Nxi.'*Nxi));
    K0e=K0e+H(iXi)*detJac*( B1.'*C*B1);
    K1e=K1e+H(iXi)*detJac*( B1.'*C*B2-B2.'*C*B1);
    K2e=K2e+H(iXi)*detJac*(-B2.'*C*B2);
end
switch dofs
    case 1
        dof = [3 6 9 12 15 18 21 24];
    case 2
        dof=[1 3 4 6 7 9 10 12 13 15 16 18 19 21 22 24];  % choose the dofs that are related, e.g., if it is 3 dofs, then include all of them (12 elements) 
    case 3
        dof=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
end
Me =Me(dof,dof);
K0e=K0e(dof,dof);
K1e=K1e(dof,dof);
K2e=K2e(dof,dof);
K3e=K3e(dof,dof);
K4e=K4e(dof,dof);
C0e=C0e(dof,dof);