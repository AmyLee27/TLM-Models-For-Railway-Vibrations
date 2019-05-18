function [Ks]=stiffvk(knum,b,uzk)

%%%% Code for the interface - 1-interaction on the middle point;
%%%% 2-interaction on the medium displacement
                              
inte = 1;    % interaction type (1=middle, 2=average) NB: 2 is more accurate for ballasted model, 1 for slab one

% Preparation of data;

k1=[knum];
k2=k1;

% Prea-allocation of variables;
uzk2=sparse(length(k2));
hs=sparse(length(k1));
Ks=sparse(length(k1));


% Computation of the dynamic stiffness

for j=1:length(k1);
    for j2=1:length(k2);
        k=k2(j2);
        if inte==1;
             if k==0;
                uzk2(j2)=uzk(j,j2);
             else
                uzk2(j2)=uzk(j,j2)*sin(k*b)/(k*b);
             end
        else
           if k==0;
                uzk2(j2)=uzk(j,j2);
             else
                uzk2(j2)=uzk(j,j2)*(sin(k*b))^2/(k*b)^2;
               
           end 
        end
    end
    hs(j)=trapz(k2,uzk2)/(2*pi);
    Ks(j)=1/hs(j);
end

