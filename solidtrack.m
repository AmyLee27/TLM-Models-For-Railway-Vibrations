function [b_node,b_element]=solidtrack(l1,l2,n1,h,n2,width_t,z0,l_rp,code)

%   SOLID mesh.
%   [b_node,b_element]=ballast(l1,l2,n1,h,n2,width_t)
%   returns the ballast mesh
%
%   l1  lower length
%   l2  upper length
%   n1  horizontal divisions
%   h   height
%   n2  vertical divisions
%   width_t track gauge/2
%   z0 vertical coordinate at the bottom
%   code 1 if sleeper

nt=(n1+1)*(n2+1);        %total number of nodes
b_node=zeros(nt,3);
b_node(:,1)=1:1:nt;
ne=n1*n2;                %total number of elements
b_element=zeros(ne,5);
b_element(:,1)=1:1:ne;

for inode=1:n1+1
    b_node(inode,2)=-l1/2+(inode-1)*l1/n1;
    b_node((n1+1)*n2+inode,2)=-l2/2+(inode-1)*l2/n1;
    b_node((n1+1)*n2+inode,3)=h;
    alpha(inode)=atan(h/(b_node((n1+1)*n2+inode,2)-b_node(inode,2)));
    if n2>1
        for in2=2:n2
                b_node((n1+1)*(in2-1)+inode,3)=(h/n2)*(in2-1);
                b_node((n1+1)*(in2-1)+inode,2)=b_node(inode,2)+(b_node((n1+1)*(in2-1)+inode,3))/tan(alpha(inode));
        end
    end
end

if code==1
b_node2=b_node(((n1+1)*n2+1):((n1+1)*n2+n1+1),:,:);
xb_node2=b_node2(:,2);
xb_node2=(xb_node2)+(1*width_t);
[minxb,indmin]=min(abs(xb_node2));
ind1=indmin;
ind2=length(xb_node2)+1-ind1;
nind1=b_node2(ind1,1);
nind2=b_node2(ind2,1);
b_node(nind1,2)=-(1*width_t);
b_node(nind2,2)=(1*width_t);

b_node(nind1-1,2)=-(1*width_t+l_rp/2);
b_node(nind1+1,2)=-(1*width_t-l_rp/2);
b_node(nind2-1,2)=(1*width_t-l_rp/2);
b_node(nind2+1,2)=(1*width_t+l_rp/2);
end
%
%
for in2=1:n2
    for inode=1:n1
        b_element((in2-1)*n1+inode,2)=(in2-1)*(n1+1)+inode;
        b_element((in2-1)*n1+inode,3)=(in2-1)*(n1+1)+1+inode;
        b_element((in2-1)*n1+inode,4)=in2*(n1+1)+1+inode;
        b_element((in2-1)*n1+inode,5)= in2*(n1+1)+inode;
    end
end

b_node(:,3)=b_node(:,3)+z0;