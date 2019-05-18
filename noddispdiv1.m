% FEM 13D - Finite element code for infinite structures

%Construction of the global matrix and computation of the
%%%%%%%%%%%%%%%%%%%%%%%%%displacements %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ax,ay,az]=noddispdiv1(k1,k2,w,nodes1,cinf,properinf1,nodinf1,K0f,K1f,K2f,K3f,K4f,K5f,M,Kr,F)
%Initialization


%%%%%%%Computation of the stiffmatrix%%%%%%%%%%%%%%%%%%%%%%
if w>0,
Kt=K0f-1i*k1*K1f-1i*k2*K2f+k1^2*K3f+k2^2*K4f+k1*k2*K5f-w^2*M+Kr;
else
if w==0,
Kt=real(K0f)-1i*k1*real(K1f)-1i*k2*real(K2f)+k1^2*real(K3f)+k2^2*real(K4f)+k1*k2*real(K5f)+Kr;
else
Kt=conj(K0f)-1i*k1*conj(K1f)-1i*k2*conj(K2f)+k1^2*conj(K3f)+k2^2*conj(K4f)+k1*k2*conj(K5f)-w^2*M+Kr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%Infinite%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%elements%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  if cinf==1,
          dof=zeros(1,3);
          Kdyn=Kdyeleinf1(properinf1,k1,k2,w);
           
          
             dof(1)=3*nodinf1-2;
             dof(2)=dof(1)+1;
             dof(3)=dof(1)+2;
        
  
     Kt(dof,dof)=Kt(dof,dof)+Kdyn;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%Solving the system%%%%%%%%%%%%%%%
    a=Kt\F;
    
      
    a=a(1:3*size(nodes1,1),1); %%%%Displacements vector%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%Computing stresses and strains%%%%%%%%%%%%%%%%%%%%%%
    [ax,ay,az]=displac(a); %Organizes the vector of displacements in 3 vectors x, y, and z
    

   
   
    
  