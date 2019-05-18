% ITM-FEM 2.5D - Finite element code for infinite structures

%Code developed by Pedro Alves Costa
% 04/07/2008
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Computes the matrix Kdyn %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kdyn=[Tq]*[G]^-1*[T]


function [Kdyn,G1]=KdynMatrix_PAC(uxxk,uyxk,uzxk,uxyk,uyyk,uzyk,uxzk,uyzk,uzzk,k2,jk1,nitmnode,zeroNod,feNod)



%%%%Computes de G matrix%%%%%%
[G]=Gmatrix_PAC(uxxk,uyxk,uzxk,uxyk,uyyk,uzyk,uxzk,uyzk,uzzk,k2,jk1,nitmnode,zeroNod,feNod);


%%%%Computes T and Tq matrix%%%%%
[T,Tq]=MatrixTTq_PAC(feNod,zeroNod,nitmnode);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kdyn=Tq*inv(G)*T;
G1=G;
end

