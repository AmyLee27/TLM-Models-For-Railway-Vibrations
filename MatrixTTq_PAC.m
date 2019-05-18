% ITM-FEM 2.5D - Finite element code for infinite structures

%Code developed by Pedro Alves Costa
% 12/06/2008
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Compoutes the matrices Tq and T %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kdyn=[Tq]*[G]^-1*[T]


function [T,Tq]=MatrixTTq_PAC(feNod,zeroNod,nitmnode)

%%%%% Reading data

%load nodes.dat; % list of the number and coordinates of the nodes of FEM mesh 
%nodes(:,1)=[];
%load nfemnode.dat; % list of the number of the nodes of the FEM mesh in contact with ITM
%load nitmnode.dat; % list of the number of the collocation nodes in ITM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%Preallocation

T=sparse(3*length(nitmnode),3*length(zeroNod));
Tq=sparse(3*length(zeroNod),3*length(nitmnode));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555

N=[0,0,0,1,0,0,0,0,0;
   0,0,0,0,1,0,0,0,0;
   0,0,0,0,0,1,0,0,0]; % Matrix N in the collocation nodes

Nint=[1/6,0,0;
      0,1/6,0;
      0,0,1/6;
      2/3,0,0;
      0,2/3,0;
      0,0,2/3;
      1/6,0,0;
      0,1/6,0;
      0,0,1/6];  % Matrix N' integrated

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for r=1:length(nitmnode)
    T(3*r-2:3*r,6*r-5:6*r+3)=N;   % Matrix T
    Tq(6*r-5:6*r+3,3*r-2:3*r)=Nint*abs(feNod(zeroNod(2*r-1),2)-feNod(zeroNod(2*r+1),2)); % Matrix Tq 
end
      
      
