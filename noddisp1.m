% FEM 13D - Finite element code for infinite structures

%Code developed by Pedro Alves Costa
% 11/02/2008 - 12/02/2008
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Computes the nodal displacements %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a,ax,ay,az]=noddisp1(k1,k2,w,nodes1,elements1,properties1,suports1,loads1,cinf,properinf1,nodinf1)
%Initialization



 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CODE INITIALIZATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%Preallocation
 Kt=sparse(3*size(nodes1,1),3*size(nodes1,1));

 
 %Assembly Procedure
 
 %%%%%%%%%%%%%%%Finite%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%elements%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 pcont=0;  %%%%% contador para ter em conta a existência de elementos iguais
  for j=1:size(elements1,1),
      pconel=properties1(j,5); % parâmetro identificador do elemento
      if pconel~=pcont,
     % Local Dynamic stifness matrix
      Kdyn=Kdynele1(nodes1(elements1(j,:),:),properties1(j,:),k1,k2,w); %matriz de 24 por 24-3*8x3*8
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
      Kt(dof,dof)=Kt(dof,dof)+Kdyn;
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
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
  %%%%Construction of the load vector - "Loads+imposed displacements on the supports"%%%%%%
   
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
    
    
    
    
    
    %%%% Construction of the matrix "K-solve-system"
    Ks=Kt;
    for j=1:length(dofi),
        Ks(dofi(j),3*size(nodes1,1)+j)=1;
        Ks(3*size(nodes1,1)+j,dofi(j))=1;
    end
    
    
    
    
    %%%%%%%%%Solving the system%%%%%%%%%%%%%%%
    a=Ks\F;
    
    a=a(1:3*size(nodes1,1),1); %%%%Displacements vector%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%Computing stresses and strains%%%%%%%%%%%%%%%%%%%%%%
    [ax,ay,az]=displac(a); %Organizes the vector of displacements in 3 vectors x, y, and z
    
    %[strt,strest]=Strain(ax,ay,az,elements,nodes,properties,k,w);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
    
  