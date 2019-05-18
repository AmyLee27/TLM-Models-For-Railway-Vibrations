function [Me,K0e,K1e,K2e,K3e,K4e,C0e]=link2(kzz,czz)

K0e=[kzz  -kzz    
    -kzz   kzz];

C0e=[czz  -czz    
    -czz   czz];

Me=zeros(2,2);
K1e=zeros(2,2);
K2e=zeros(2,2);
K3e=zeros(2,2);
K4e=zeros(2,2);

