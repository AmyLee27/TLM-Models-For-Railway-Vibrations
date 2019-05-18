function [b_node,b_element]=solidtrack_8nodes_new(l1,l2,n1,h,n2,width_t,z0,l_rp,code)

%   SOLID mesh.
%   [b_node,b_element]=ballast(l1,l2,n1,h,n2,width_t)
%   returns the ballast mesh
%
%   l1  lower length
%   l2  upper length
%   n1  horizontal divisions - number of horizontal elements
%   h   height
%   n2  vertical divisions - number of vertical elements
%   width_t track gauge/2
%   z0 vertical coordinate at the bottom
%   code 1 if sleeper / code 2 if it is 7-node element


if code ==0 || code ==1  
    nt=(2*n1+1)*(n2+1)+n2*(n1+1);        %total number of nodes
    ne=n1*n2;                %total number of elements  
elseif code ==2   
    nt=(2*n1+1)*n2 + (n2+1)*(n1+1);        %total number of nodes
    ne=n1*n2;                %total number of elements  
end

b_node=zeros(nt,3);
b_node(:,1)=1:1:nt;
b_node3=zeros((2*n1+1)*(n2+1)+n2*(n1+1),3);
b_node3(:,1)=1:1:(2*n1+1)*(n2+1)+n2*(n1+1);
b_element=zeros(ne,9);
b_element(:,1)=1:1:ne;
 
if code ==0 || code ==1 
    for inode=1:(2*n1+1)
        b_node(inode,2)=-l1/2+(inode-1)*l1/(2*n1);
        b_node((2*n1+1)*n2+n2*(n1+1)+inode,2)=-l2/2+(inode-1)*l2/(2*n1);
        b_node((2*n1+1)*n2+n2*(n1+1)+inode,3)=h;
        alpha(inode)=atan(h/(b_node((2*n1+1)*n2+n2*(n1+1)+inode,2)-b_node(inode,2)));

        if n2>1
            for in2=2:n2
                b_node((2*n1+1)*(in2-1)+(in2-1)*(n1+1)+inode,3)=(h/n2)*(in2-1);
                b_node((2*n1+1)*(in2-1)+(in2-1)*(n1+1)+inode,2)=b_node(inode,2)+(b_node((2*n1+1)*(in2-1)+(in2-1)*(n1+1)+inode,3))/tan(alpha(inode));
            end
        end
    end

    for in3=1:n2
        b_node(((2*n1+1)*in3+(n1+1)*(in3-1)+1):((2*n1+1)*in3+(n1+1)*(in3-1)+n1+1),3)=h/(2*n2)*(2*in3-1);
        b_node(((2*n1+1)*in3+(n1+1)*(in3-1)+1):((2*n1+1)*in3+(n1+1)*(in3-1)+n1+1),2)=0.50*(b_node((1+(3*n1+2)*(in3-1)):2:(2*n1+1+(3*n1+2)*(in3-1)),2)+b_node((1+(3*n1+2)*in3):2:(2*n1+1+(3*n1+2)*in3),2));
    end

elseif code ==2
    
    for inode=1:(2*n1+1)
        b_node3(inode,2)=-l1/2+(inode-1)*l1/(2*n1);
        b_node3((2*n1+1)*n2+n2*(n1+1)+inode,2)=-l2/2+(inode-1)*l2/(2*n1);
        b_node3((2*n1+1)*n2+n2*(n1+1)+inode,3)=h;
        alpha(inode)=atan(h/(b_node3((2*n1+1)*n2+n2*(n1+1)+inode,2)-b_node3(inode,2)));

        if n2>1
            for in2=2:n2
                b_node3((2*n1+1)*(in2-1)+(in2-1)*(n1+1)+inode,3)=(h/n2)*(in2-1);
                b_node3((2*n1+1)*(in2-1)+(in2-1)*(n1+1)+inode,2)=b_node3(inode,2)+(b_node3((2*n1+1)*(in2-1)+(in2-1)*(n1+1)+inode,3))/tan(alpha(inode));
            end
        end
    end

    for in3=1:n2
        b_node3(((2*n1+1)*in3+(n1+1)*(in3-1)+1):((2*n1+1)*in3+(n1+1)*(in3-1)+n1+1),3)=h/(2*n2)*(2*in3-1);
        b_node3(((2*n1+1)*in3+(n1+1)*(in3-1)+1):((2*n1+1)*in3+(n1+1)*(in3-1)+n1+1),2)=0.50*(b_node3((1+(3*n1+2)*(n2-1)):2:(2*n1+1+(3*n1+2)*(n2-1)),2)+b_node3((1+(3*n1+2)*n2):2:(2*n1+1+(3*n1+2)*n2),2));
    end
    
    b_node1(:,:)=b_node3(1:2:(2*n1+1),:);
    b_node1(:,1) = 1:1:(n1+1);
    b_node = [b_node1 ;b_node3(2*n1+2:end,:)];
    b_node(:,1) = 1:1:nt;
end

if code==1
b_node2=b_node(((2*n1+1)*n2+1+n2*(n1+1)):((2*n1+1)*(n2+1)+n2*(n1+1)),:,:);
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
if code == 1 || code == 0
    for in2=1:n2
        for inode=1:n1

            b_element((in2-1)*n1+inode,2)=(in2-1)*(2*n1+1)+(in2-1)*(n1+1)+2*inode-1;%1
            b_element((in2-1)*n1+inode,3)=(in2-1)*(2*n1+1)+(in2-1)*(n1+1)+2*inode+1;%3
            b_element((in2-1)*n1+inode,4)=in2*(3*n1+2)+2*inode+1;%5
            b_element((in2-1)*n1+inode,5)=in2*(3*n1+2)+2*inode-1;%7
            b_element((in2-1)*n1+inode,6)=(in2-1)*(2*n1+1)+(in2-1)*(n1+1)+2*inode;%2
            b_element((in2-1)*n1+inode,7)=in2*(2*n1+1)+(in2-1)*(n1+1)+inode+1;%4
            b_element((in2-1)*n1+inode,8)=in2*(3*n1+2)+2*inode;%6
            b_element((in2-1)*n1+inode,9)= in2*(2*n1+1)+(in2-1)*(n1+1)+inode;%8                

%             b_element((in2-1)*n1+inode,2)=(in2-1)*(2*n1+1)+(in2-1)*(n1+1)+2*inode-1;%1
%             b_element((in2-1)*n1+inode,3)=(in2-1)*(2*n1+1)+(in2-1)*(n1+1)+2*inode;%2
%             b_element((in2-1)*n1+inode,4)=(in2-1)*(2*n1+1)+(in2-1)*(n1+1)+2*inode+1;%3
%             b_element((in2-1)*n1+inode,5)=in2*(2*n1+1)+(in2-1)*(n1+1)+inode+1;%4
%             b_element((in2-1)*n1+inode,6)=in2*(3*n1+2)+2*inode+1;%5
%             b_element((in2-1)*n1+inode,7)=in2*(3*n1+2)+2*inode;%6
%             b_element((in2-1)*n1+inode,8)=in2*(3*n1+2)+2*inode-1;%7
%             b_element((in2-1)*n1+inode,9)= in2*(2*n1+1)+(in2-1)*(n1+1)+inode;%8 
                       
        end
    end
elseif code == 2
      
    for inode=1:n1

        b_element(inode,2)=inode;%1
        b_element(inode,3)=inode+1;
        b_element(inode,4)=n1+1+inode+1;
        b_element(inode,5)=2*(n1+1)+2*inode+1;
        b_element(inode,6)=2*(n1+1)+2*inode;
        b_element(inode,7)=2*(n1+1)+2*inode-1;
        b_element(inode,8)=n1+1+inode;
        b_element(inode,9)=0;               

    end
    
    for in2=2:n2
        for inode=1:n1

            b_element((in2-1)*n1+inode,2)=(n1+1)+(n1+1)*(in2-1)+(2*n1+1)*(in2-2)+2*inode-1;%1
            b_element((in2-1)*n1+inode,3)=(n1+1)+(n1+1)*(in2-1)+(2*n1+1)*(in2-2)+2*inode;%2
            b_element((in2-1)*n1+inode,4)=(n1+1)+(n1+1)*(in2-1)+(2*n1+1)*(in2-2)+2*inode+1;%3
            b_element((in2-1)*n1+inode,5)=(n1+1)+(n1+1)*(in2-1)+(2*n1+1)*(in2-1)+inode+1;%4
            b_element((in2-1)*n1+inode,6)=(n1+1)+(n1+1)*in2+(2*n1+1)*(in2-1)+2*inode+1;%5
            b_element((in2-1)*n1+inode,7)=(n1+1)+(n1+1)*in2+(2*n1+1)*(in2-1)+2*inode;%6
            b_element((in2-1)*n1+inode,8)=(n1+1)+(n1+1)*in2+(2*n1+1)*(in2-1)+2*inode-1;%7
            b_element((in2-1)*n1+inode,9)=(n1+1)+(n1+1)*(in2-1)+(2*n1+1)*(in2-1)+inode;%8                              

        end
    end
    
    
end
b_node(:,3)=b_node(:,3)+z0;