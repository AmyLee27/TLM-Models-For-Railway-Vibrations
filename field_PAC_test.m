% calculate free field response

function [uf3,s]=field_PAC_test(input_x,input_y,input_z,k2,k1,nitmnode,zeroNod,feNod,y0,traction_x,traction_y,traction_z)

%%%%%Atenção jk1=indicador da posição de k1 em estudo
npoi=length(nitmnode);    
ny=length(y0);
nd=length(input_x(1,1,:));
nk1=length(k1);

uf1 = 0;
uf2 = zeros(ny,nk1);
uf3 = zeros(nk1,ny,nd);

for l = 1:nd
    for jk1=1:nk1
        for ir=1:ny      
            for ij=1:npoi      %%%%ATENÇÃO - A FORÇA ESTÁ APLICADA NO NÓ IJ
                le=abs(feNod(zeroNod(2*ij-1),2)-feNod(zeroNod(2*ij+1),2));
                R_x=2*traction_x(ij,jk1).*(sin(k2*le/2)./(k2)).*exp(1i*k2*(y0(ir)-nitmnode(ij)));
                R_x((length(k2)+1)/2)=le.*traction_x(ij,jk1);
                
                R_y=2*traction_y(ij,jk1).*(sin(k2*le/2)./(k2)).*exp(1i*k2*(y0(ir)-nitmnode(ij)));
                R_y((length(k2)+1)/2)=le.*traction_y(ij,jk1);
                
                R_z=2*traction_z(ij,jk1).*(sin(k2*le/2)./(k2)).*exp(1i*k2*(y0(ir)-nitmnode(ij)));
                R_z((length(k2)+1)/2)=le.*traction_z(ij,jk1);

                u=trapz(k2,(R_x.*input_x(jk1,:,l)+R_y.*input_y(jk1,:,l)+R_z.*input_z(jk1,:,l)))/(2*pi);
                uf1 = uf1 + u;                
            end
    
            uf2(ir,jk1) = uf1;
            uf1 = 0;
        end
    end
    
    for ir=1:ny
        [s,transform]=inverse(k1,uf2(ir,:)); 
        uf3(:,ir,l)=transform;
    end
end




end