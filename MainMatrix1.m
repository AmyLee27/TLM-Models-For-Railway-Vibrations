% FEM 13D - Finite element code for infinite structures
%%%Computation of the main matrices%%%%%%%%%%%%%%
%Code developed by Pedro Alves Costa
%%%20/10/2008
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Computes the nodal displacements %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F]=MainMatrix1(nodes1,elements1,properties1,suports1,loads1)
%Initialization



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CODE INITIALIZATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%Preallocation
 K0f=sparse(3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)),3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)));
 K1f=sparse(3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)),3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)));
 K2f=sparse(3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)),3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)));
 K3f=sparse(3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)),3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)));
 K4f=sparse(3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)),3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)));
 K5f=sparse(3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)),3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)));
 Mf=sparse(3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)),3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)));
 Kr=sparse(3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)),3*size(nodes1,1)+sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4)));

 
 %Assembly Procedure
 
 %%%%%%%%%%%%%%%Finite%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%elements%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 pcont=0;  %%%%% contador para ter em conta a existência de elementos iguais
 %  counter to take account of the existence of the same elements
  for j=1:size(elements1,1),
      pconel=properties1(j,5); % parâmetro identificador do elemento
                               % parameter identifier element
      if pconel~=pcont,
     % Local Dynamic stifness matrix
     
     [K0,K1,K2,K3,K4,K5,M]=KMelematrixdiv1(nodes1(elements1(j,:),:),properties1(j,:));
      dof=zeros(1,3*length(elements1(1,:)));
      pcont=pconel;
      end
      
      % Allocation to the correct position
     
     % Identification of the numeration of the degrees of freedom
        ele=elements1(j,:);
        for r=1:length(ele),
             dof(3*r-2)=3*ele(r)-2;
             dof(3*r-1)=dof(3*r-2)+1;
             dof(3*r)=dof(3*r-2)+2;
        end
      
      %Assembly
      K0f(dof,dof)=K0f(dof,dof)+K0;
      K1f(dof,dof)=K1f(dof,dof)+K1;
      K2f(dof,dof)=K2f(dof,dof)+K2;
      K3f(dof,dof)=K3f(dof,dof)+K3;
      K4f(dof,dof)=K4f(dof,dof)+K4;
      K5f(dof,dof)=K5f(dof,dof)+K5;
      Mf(dof,dof)=Mf(dof,dof)+M;
      
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
%%%%%%%%Construction of the vector loads%%%%%%%%%%%%%%%%%
    %%%%Preallocation
   dofimp=sum(suports1(:,2))+sum(suports1(:,3))+sum(suports1(:,4));
   F=zeros(3*size(nodes1,1)+dofimp,1);
   dofi=sparse(dofimp,1); 
   
    %%%% Run the vector loads%%%%%%%
    for j=1:size(loads1,1),
            if loads1(j,2)~=0,                    %%%%load with the x direction%%%%%%%%
                F(3*loads1(j,1)-2,1)=loads1(j,2);  
            end
            if loads1(j,3)~=0,
                F(3*loads1(j,1)-1,1)=loads1(j,3);          %%%%load with the y direction%%%%%%%%
            end
            if loads1(j,4)~=0,
                F(3*loads1(j,1),1)=loads1(j,4);          %%%%load with the z direction%%%%%%%%
            end
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 %%%%Construction of the Kr matrix%%%%%%%%%%%
     %%%% Run the vector suports%%%%%
     ir=0;
    for j=1:size(suports1,1),
            if suports1(j,2)==1,
                ir=ir+1;
                dofi(ir)=3*suports1(j,1)-2;    %%%% support with the x direction %%%
            end
            if suports1(j,3)==1,
                ir=ir+1;
                dofi(ir)=3*suports1(j,1)-1;    %%%% support with the y direction %%%
            end
            if suports1(j,4)==1,
                ir=ir+1;
                dofi(ir)=3*suports1(j,1);      %%%% support with the z direction %%%%%
            end
    end
    dofi=sort(dofi);

    for j=1:length(dofi),
        Kr(dofi(j),3*size(nodes1,1)+j)=1;
        Kr(3*size(nodes1,1)+j,dofi(j))=1;
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
   