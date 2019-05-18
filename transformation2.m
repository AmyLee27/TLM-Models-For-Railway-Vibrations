function [post_variable,s] = transformation2(variable,nFreq,py,yy,base_nNod,sigma_mat,node_space,feNod)

uf=0;
R_scale=zeros(length(yy),length(py));
post_transform=zeros(length(py),length(yy),length(variable(1,1,:,1)));
post_variable=zeros(length(py),length(yy),length(variable(1,1,:,1)),nFreq);
for iFreq=1:nFreq
        variable1=variable(:,:,:,iFreq);
        for l=1:length(variable(1,1,:,1))
            for j=1:length(py)
                variable2=variable1(j,:,l);

                for i=1:length(yy)
                    for iNod=1:base_nNod
                        T=[0 0 1]*sigma_mat(:,:,iNod,iFreq)*[0 0 1]';
                        P_scale=2*variable2.*T.*(sin(py*node_space)./py); 
                        ii=find(py==0);
                        P_scale(ii)=2*variable2(ii).*T;
                        uf=uf+trapz(py,P_scale.*exp(1i*py*(yy(i)-feNod(iNod,2))))/(2*pi);
                    end

                    R_scale(i,j) = uf;
                end
                for i = 1:length(yy)
                    [s,transform]=inverse(py,R_scale(i,:)); 
                    post_transform(:,i,l)=transform;
                end
            
            end
        end
        post_variable(:,:,:,iFreq) = post_transform;
end