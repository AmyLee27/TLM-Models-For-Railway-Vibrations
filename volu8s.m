function [epsilon,sigma]=volu8s(x,u,du,E,nu)

lambda=nu*E/(1+nu)/(1-2*nu);
mu=E/2/(1+nu);

C=[lambda+2*mu lambda      lambda      0  0  0
    lambda      lambda+2*mu lambda      0  0  0
    lambda      lambda      lambda+2*mu 0  0  0
    0           0           0           mu 0  0
    0           0           0           0  mu 0
    0           0           0           0  0  mu];

xi=[-1 -1;
     0 -1;
     1  -1;
     1  0;
     1  1;
     0  1;
     -1 1;
    -1 0];

nXi=8;

N=zeros(nXi,8);
N(:,1)=0.25*(1.0-xi(:,1)).*(xi(:,2)-1.0).*(xi(:,1)+xi(:,2)+1.0);
N(:,2)=0.5*(1.0-xi(:,2)).*(1.0-xi(:,1).^2);
N(:,3)=0.25*(1.0+xi(:,1)).*(xi(:,2)-1.0).*(-xi(:,1)+xi(:,2)+1.0);
N(:,4)=0.5*(1.0+xi(:,1)).*(1.0-xi(:,2).^2);
N(:,5)=0.25*(1.0+xi(:,1)).*(1.0+xi(:,2)).*(xi(:,1)+xi(:,2)-1.0);
N(:,6)=0.5*(1.0+xi(:,2)).*(1.0-xi(:,1).^2);
N(:,7)=-0.25*(1.0-xi(:,1)).*(1.0+xi(:,2)).*(xi(:,1)-xi(:,2)+1.0);
N(:,8)=0.5*(1.0-xi(:,1)).*(1.0-xi(:,2).^2);

dN=zeros(nXi,16);
dN(:,1)=0.25*(1.0-xi(:,2)).*(2*xi(:,1)+xi(:,2));
dN(:,9)=0.25*(1.0-xi(:,1)).*(2*xi(:,2)+xi(:,1));
dN(:,2)=-1.0*xi(:,1).*(1.0-xi(:,2));
dN(:,10)=-0.5*(1.0-xi(:,1).^2);
dN(:,3)= -0.25*(1.0-xi(:,2)).*(xi(:,2)-2*xi(:,1));
dN(:,11)=0.25*(1.0+xi(:,1)).*(2*xi(:,2)-xi(:,1));
dN(:,4)=0.5*(1.0-xi(:,2).^2);
dN(:,12)=-1.0*xi(:,2).*(1.0+xi(:,1));
dN(:,5)= 0.25*(1.0+xi(:,2)).*(2*xi(:,1)+xi(:,2));
dN(:,13)=0.25*(1.0+xi(:,1)).*(2*xi(:,2)+xi(:,1));
dN(:,6)=-1.0*xi(:,1).*(1.0+xi(:,2));
dN(:,14)= 0.5*(1.0-xi(:,1).^2);
dN(:,7)=-0.25*(1.0+xi(:,2)).*(xi(:,2)-2*xi(:,1));
dN(:,15)=0.25*(1.0-xi(:,1)).*(2*xi(:,2)-xi(:,1));
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

epsilon=zeros(6,8);
sigma=zeros(6,8);

for iXi=1:nXi
  JacUtil=[Jac(iXi,1,1),Jac(iXi,1,2);Jac(iXi,2,1),Jac(iXi,2,2)];
  DN1= JacUtil\[dN(iXi,1);dN(iXi,9)];% [dN1/dx; dN1/dz]
  DN2= JacUtil\[dN(iXi,2);dN(iXi,10)];% [dN2/dx; dN2/dz]
  DN3= JacUtil\[dN(iXi,3);dN(iXi,11)];% [dN3/dx; dN3/dz]
  DN4= JacUtil\[dN(iXi,4);dN(iXi,12)];% [dN4/dx; dN4/dz]
  DN5= JacUtil\[dN(iXi,5);dN(iXi,13)];
  DN6= JacUtil\[dN(iXi,6);dN(iXi,14)];
  DN7= JacUtil\[dN(iXi,7);dN(iXi,15)];
  DN8= JacUtil\[dN(iXi,8);dN(iXi,16)];

  
  epsilon(1,iXi) = DN1(1)*u(1)+DN2(1)*u(4)+DN3(1)*u(7)+DN4(1)*u(10)+DN5(1)*u(13)+DN6(1)*u(16)+DN7(1)*u(19)+DN8(1)*u(22); % epsilon_xx
  epsilon(2,iXi) = du(2)+du(5)+du(8)+du(11)+du(14)+du(17)+du(20)+du(23); %strain yy
  epsilon(3,iXi)=DN1(2)*u(3)+DN2(2)*u(6)+DN3(2)*u(9)+DN4(2)*u(12)+DN5(2)*u(15)+DN6(2)*u(18)+DN7(2)*u(21)+DN8(2)*u(24);  %epsilon_zz
  epsilon(4,iXi)= (DN1(1)*u(2)+DN2(1)*u(5)+DN3(1)*u(8)+DN4(1)*u(11)+DN5(1)*u(14)+DN6(1)*u(17)+DN7(1)*u(20)+DN8(1)*u(23)+du(1)+du(4)+du(7)+du(10)+du(13)+du(16)+du(19)+du(22));% strain xy
  epsilon(5,iXi)=(DN1(1)*u(3)+DN2(1)*u(6)+DN3(1)*u(9)+DN4(1)*u(12)+DN5(1)*u(15)+DN6(1)*u(18)+DN7(1)*u(21)+DN8(1)*u(24)+DN1(2)*u(1)+DN2(2)*u(4)+DN3(2)*u(7)+DN4(2)*u(10)+DN5(2)*u(13)+DN6(2)*u(16)+DN7(2)*u(19)+DN8(2)*u(22));  %epsilon_xz
  epsilon(6,iXi)=(du(3)+du(6)+du(9)+du(12)+du(15)+du(18)+du(21)+du(24)+DN1(2)*u(2)+DN2(2)*u(5)+DN3(2)*u(8)+DN4(2)*u(11)+DN5(2)*u(14)+DN6(2)*u(17)+DN7(2)*u(20)+DN8(2)*u(22));  % strain yz                                        %epsilon_zy 
  
%   sigma(1,iXi)=epsilon(1,iXi)*(lambda+2*mu);                       %sigma_zz
%   sigma(2,iXi)=epsilon(1,iXi)*mu;                                  %sigma_xz
%   sigma(3,iXi)=epsilon(1,iXi)*mu;                                  %sigma_zy 

  sigma(:,iXi) = C*epsilon(:,iXi);

end

end