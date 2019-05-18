function [Ks]=stiffnbvk(k1max,k1inc,b,uzk)

%%%% npoi - number of interaction points
                              
npoi = 6;    % how many elements to discretise (e.g.1=same as prev case(middle),
             % if large then require large wavenumbers, good number is around=3/0.5=6)
             
%call vibfemgreen

% Preparation of data;


k1=-k1max:k1inc:k1max;
k2=k1;
Poin=1;

% Pre-allocation of variables;

H=sparse(length(npoi),length(npoi));
Ks=sparse(length(k1));
Ksr=sparse(length(k1),npoi);


%%% Discretization of the interface

dy=2*b/npoi; %%%% dimension between two consecutive intercation points

 r=1:npoi;
 yp=-b+dy*(r.*2-1)/2; %%% vector with the interaction points
    un=sparse(npoi,1);
    un(:,1)=1;



% Computation of the dynamic stiffness

for j=1:length(k1);
    
    
    for ir=1:npoi,
        
        for ij=1:npoi,
            
                  R=(sin(k2*dy/2)./(k2*dy/2)).*exp(i*k2*(yp(ir)-yp(ij)));
                  R((length(k2)+1)/2)=1;
                                           
                  uzk2=R.*uzk(j,:,Poin);
                      
             
   
            H(ir,ij)=trapz(k2,uzk2)/(2*pi);
        end
    end
    Xs=inv(H)*un;
   for ri=1:npoi,
       Ksr(j,ri)=Xs(ri,1);
       Ks(j)=sum(Ksr(j,:));
   end
end

