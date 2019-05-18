function [epsilon,sigma]=volu4s(x,u,du,E,nu)

lambda=nu*E/(1+nu)/(1-2*nu);
mu=E/2/(1+nu);

C=[lambda+2*mu lambda      lambda      0  0  0
    lambda      lambda+2*mu lambda      0  0  0
    lambda      lambda      lambda+2*mu 0  0  0
    0           0           0           mu 0  0
    0           0           0           0  mu 0
    0           0           0           0  0  mu];

xi=[-1 -1;
     1 -1;
     1  1;
    -1 1];
nXi=4;

N=zeros(nXi,4);
N(:,1)=0.25*(1.0-xi(:,1)).*(1.0-xi(:,2));                     
N(:,2)=0.25*(1.0+xi(:,1)).*(1.0-xi(:,2));                    
N(:,3)=0.25*(1.0+xi(:,1)).*(1.0+xi(:,2));                    
N(:,4)=0.25*(1.0-xi(:,1)).*(1.0+xi(:,2));                     

dN=zeros(nXi,8);
dN(:,1)=-0.25*(1.0-xi(:,2));                                 
dN(:,2)= 0.25*(1.0-xi(:,2));                                 
dN(:,3)= 0.25*(1.0+xi(:,2));                                 
dN(:,4)=-0.25*(1.0+xi(:,2));                                 
dN(:,5)=-0.25*(1.0-xi(:,1));                                
dN(:,6)=-0.25*(1.0+xi(:,1));                                 
dN(:,7)= 0.25*(1.0+xi(:,1));                                  
dN(:,8)= 0.25*(1.0-xi(:,1));   

Jac=zeros(nXi,2,2);
for iXi=1:nXi
	for iNod=1:4
	  Jac(iXi,1,1)=Jac(iXi,1,1)+dN(iXi,iNod)*x(iNod,1);
	  Jac(iXi,1,2)=Jac(iXi,1,2)+dN(iXi,iNod)*x(iNod,2);
	  Jac(iXi,2,1)=Jac(iXi,2,1)+dN(iXi,4+iNod)*x(iNod,1);
	  Jac(iXi,2,2)=Jac(iXi,2,2)+dN(iXi,4+iNod)*x(iNod,2);
  end
end

epsilon=zeros(6,4);
sigma=zeros(6,4);

for iXi=1:nXi
  JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2);Jac(iXi,2,1),Jac(iXi,2,2)];
  DN1= JacUtil\[dN(iXi,1);dN(iXi,5)]; % [dN1/dx; dN1/dz]
  DN2= JacUtil\[dN(iXi,2);dN(iXi,6)];% [dN2/dx; dN2/dz]
  DN3= JacUtil\[dN(iXi,3);dN(iXi,7)]; % [dN3/dx; dN3/dz]
  DN4= JacUtil\[dN(iXi,4);dN(iXi,8)]; % [dN4/dx; dN4/dz]
  
  epsilon(1,iXi) = DN1(1)*u(1)+DN2(1)*u(4)+DN3(1)*u(7)+DN4(1)*u(10); % epsilon_xx
  epsilon(2,iXi) = du(2)+du(5)+du(8)+du(11); %strain yy
  epsilon(3,iXi)=DN1(2)*u(3)+DN2(2)*u(6)+DN3(2)*u(9)+DN4(2)*u(12);  %epsilon_zz
  epsilon(4,iXi)= DN1(1)*u(2)+DN2(1)*u(5)+DN3(1)*u(8)+DN4(1)*u(11)+du(1)+du(4)+du(7)+du(10);% strain xy
  epsilon(5,iXi)=DN1(1)*u(3)+DN2(1)*u(6)+DN3(1)*u(9)+DN4(1)*u(12)+DN1(2)*u(1)+DN2(2)*u(4)+DN3(2)*u(7)+DN4(2)*u(10);  %epsilon_xz
  epsilon(6,iXi)=du(3)+du(6)+du(9)+du(12)+DN1(2)*u(2)+DN1(2)*u(5)+DN1(2)*u(8)+DN1(2)*u(11);  % strain yz                                        %epsilon_zy 
  
%   sigma(1,iXi)=epsilon(1,iXi)*(lambda+2*mu);                       %sigma_zz
%   sigma(2,iXi)=epsilon(1,iXi)*mu;                                  %sigma_xz
%   sigma(3,iXi)=epsilon(1,iXi)*mu;                                  %sigma_zy 

  sigma(:,iXi) = C*epsilon(:,iXi);

end

end