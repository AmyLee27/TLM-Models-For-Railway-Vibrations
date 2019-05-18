% ITM-FEM 2.5D - Finite element code for infinite structures

%Code developed by Pedro Alves Costa
% 12/06/2008
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Computes the matrix G %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kdyn=[Tq]*[G]^-1*[T]


function [G]=Gmatrix_PAC(uxxk,uyxk,uzxk,uxyk,uyyk,uzyk,uxzk,uyzk,uzzk,k2,jk1,nitmnode,zeroNod,feNod)

%%%%%Atenção jk1=indicador da posição de k1 em estudo
npoi=length(nitmnode);    
 
for ir=1:npoi        
    for ij=1:npoi      %%%%ATENÇÃO - A FORÇA ESTÁ APLICADA NO NÓ IJ
        le=abs(feNod(zeroNod(2*ij-1),2)-feNod(zeroNod(2*ij+1),2));
        R=2*(sin(k2*le/2)./(k2)).*exp(1i*k2*(nitmnode(ir)-nitmnode(ij)));
        R((length(k2)+1)/2)=le;
                                          
        Gij(1,1)=trapz(k2,R.*uxxk(jk1,:))/(2*pi);
        Gij(1,2)=trapz(k2,R.*uxyk(jk1,:))/(2*pi);
        Gij(1,3)=trapz(k2,R.*uxzk(jk1,:))/(2*pi);
        Gij(2,1)=trapz(k2,R.*uyxk(jk1,:))/(2*pi);
        Gij(2,2)=trapz(k2,R.*uyyk(jk1,:))/(2*pi);
        Gij(2,3)=trapz(k2,R.*uyzk(jk1,:))/(2*pi);
        Gij(3,1)=trapz(k2,R.*uzxk(jk1,:))/(2*pi);
        Gij(3,2)=trapz(k2,R.*uzyk(jk1,:))/(2*pi);
        Gij(3,3)=trapz(k2,R.*uzzk(jk1,:))/(2*pi);
            
        G(3*ir-2:3*ir,3*ij-2:3*ij)=Gij;
            
    end
end

         
     %%%%%Colocar aqui a possibilidade de simetria e antisimetria    
         

     
 end