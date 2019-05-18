% calculate free field response

function [uf]=field_PAC(uzzk,k2,jk1,nitmnode,zeroNod,feNod,y0,traction)

%%%%%Atenção jk1=indicador da posição de k1 em estudo
npoi=length(nitmnode);    
ny=length(y0);
uf = 0;

for ir=1:ny      
    for ij=1:npoi      %%%%ATENÇÃO - A FORÇA ESTÁ APLICADA NO NÓ IJ
        le=abs(feNod(zeroNod(2*ij-1),2)-feNod(zeroNod(2*ij+1),2));
        R=2*traction(ij).*(sin(k2*le)./(k2)).*exp(1i*k2*(y0(ir)-nitmnode(ij)));
        R((length(k2)+1)/2)=le;

        u=trapz(k2,R.*uzzk(jk1,:))/(2*pi);
        uf = uf + u;             
    end
    
   
end
  
end