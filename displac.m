function [ax,ay,az]=displac(a)

%%%%Function to organize the displacement vector in subvectors%%%%
%
%
%%%%%%%%%%%%%%%Developed by Pedro Costa - 06/07/2007%%%%%%%%%%%%%%



%Preallocation%%%%%%%%%
ax=zeros(length(a)/3,1);
ay=zeros(length(a)/3,1);
az=zeros(length(a)/3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(a)/3,
    ax(j)=a(3*j-2);
    ay(j)=a(3*j-1);
    az(j)=a(3*j);
end