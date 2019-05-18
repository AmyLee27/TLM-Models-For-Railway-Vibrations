function [b_node,b_element]=solidtrack_7nodes(l1,l2,n1,h,n2,z0)

%   SOLID mesh.
%   [b_node,b_element]=ballast(l1,l2,n1,h,n2,width_t)
%   returns the ballast mesh
%
%   l1  lower length
%   l2  upper length
%   n1  horizontal divisions - number of horizontal elements
%   h   height
%   n2  vertical divisions - number of vertical elements (it has to be 1 for 7-node element)
%   width_t track gauge/2
%   z0 vertical coordinate at the bottom
%   code 1 if sleeper
n2=1;
nt=4*n1+3;        %total number of nodes
b_node=zeros(nt,3);
b_node(:,1)=1:1:nt;
ne=n1*n2;                %total number of elements
b_element=zeros(ne,8);
b_element(:,1)=1:1:ne;

for inode=1:(n1+1)
    b_node(inode,2)=-l1/2+(inode-1)*l1/n1;
    b_node(2*n1+2+2*inode-1,2)=-l2/2+(inode-1)*l2/n1;
    b_node(2*n1+2+2*inode-1,3)=h;
    if inode >1
        b_node(2*n1+2+2*inode-2,2) = -l2/2+(inode-1)*l2/n1 + l2/(2*n1);
        b_node(2*n1+2+2*inode-2,3) = h;
    end 
    b_node(n1+1+inode,3)=h/2;
    b_node(n1+1+inode,2)=(1/2)*(b_node(inode,2)+b_node(2*n1+2+2*inode-1,2));

end

%
%

for inode=1:n1
    
        b_element(inode,2)=inode;%1
        b_element(inode,3)=inode+1;
        b_element(inode,4)=n1+1+inode+1;
        b_element(inode,5)=2*(n1+1)+2*n1+1;
        b_element(inode,6)=2*(n1+1)+2*n1;
        b_element(inode,7)=2*(n1+1)+2*n1-1;
        b_element(inode,8)=n1+1+inode;
                
end


b_node(:,3)=b_node(:,3)+z0;