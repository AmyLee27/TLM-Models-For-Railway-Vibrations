function [US]=Reorde(u,k1)
% faz o ordenamento do vector admitindo que para k1 o resultado n é nulo
US=zeros(2*length(k1)-1,length(u(1,:)));
US(1:length(k1),:)=u;
%US(length(k1),:)=real(u(length(k1)-1,:));
for s=1:length(k1)-1;
    US(s+length(k1),:)=real(u(length(k1)-s,:))-i*imag(u(length(k1)-s,:));
end
    
    
