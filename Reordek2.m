function [US]=Reordek2(u,k2)
% faz o ordenamento do vector admitindo que para k1 o resultado é nulo
% os números de onda devem sre crescentes de -kmax ate kmax.
US=zeros(length(u(:,1)),2*length(k2)-1);
US(:,1:length(k2))=u;
%US(:,length(k2))=0;
for s=1:length(k2)-1;
    US(:,s+length(k2))=real(u(:,length(k2)-s))+1i*imag(u(:,length(k2)-s));
end
    
    
