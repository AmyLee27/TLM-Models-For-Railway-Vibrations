function [post_transform,s] = transformation(k1,k2,B,P,pre_transform,yy)

% Transform the strans/stresses calculated from the greenfunction to the
% time space domain, both in k1(-x) and k2(-y) directions

% Input
% k1 - wavenumber in x-direction
% k2 - wavenumber in y-direction
% B - Track width
% P - Contact force between the base of the track structure and soil (formed by equivalent stiffness)
% pre_transform - stress/strain calculated in green fucntion in wavenumber domain
% yy - y axis in space domain

% Output
% post_transform - stress/strain in time space domain
% s - x axis in space domain

R_scale2 = zeros(length(yy),length(k1));

for l = 1:length(pre_transform(1,1,:))

    for j=1:length(k1)
    
        P_scale=P(j)*sin(k2*B/2)./(k2*B/2);      
        ii=find(k2==0);
        P_scale(ii)=P(ii);
        R_scale=P_scale.*pre_transform(j,:,l);
                
        for m = 1:length(yy)
            R_scale2(m,j)=trapz(k2,R_scale.*exp(-1i*k2*yy(m)))/(2*pi);               
        end
        
    end
        
        for m = 1:length(yy)
            [s,transform]=inverse(k1,R_scale2(m,:)); 
            post_transform(:,m,l)=transform;
        end
       
end