function [uzk]=greenfunction(knum,cinf,v,nodes1,elements1,properties1,properinf1)

ki=knum;
kj=ki;

elestrain = 1:length(elements1(:,1));
uzk=zeros(length(ki),length(kj)); % disp infor for the top node

%Finite elements properties
no=max(nodes1(:,1));
nodes1(:,1)=[];
elements1(:,1)=[];
properties1(:,1)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%

%Infinite elements properties
if cinf==1;
    properinf1=properinf1;  % property for bottom homogeneous layer
    nodinf1=1;
else
    properinf1=0;   
    nodinf1=0;
end

% Suports and loads%%%%%
if cinf==1;
    suports1=[1,0,0,0];
else
    suports1=[1,1,1,1];
end
loads1=[length(nodes1(:,1)), 0,0,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F]=MainMatrix1(nodes1,elements1,properties1,suports1,loads1);

parfor j=1:length(ki);
    
    k1=ki(j);  % It is calculated in k1 direction, i.e., longitudinal direction (-x)
    w=-k1*v; % k(wavenumber) = 2*pi*f/v = w(angular frequency)/v;
   [uzk2]=Mnoddispdiv1(kj,k1,w,nodes1,cinf,elements1,properinf1,nodinf1,K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F,no,elestrain,properties1);

    uzk(j,:)=uzk2;

end
% uzk(length(ki),length(kj))=uzk(length(ki),length(kj)-1);
% [US]=Reordek2(uzk,kj);
% [US]=Reorde(US,ki);
% uzk=US;
% 
% na = find(ki==0);
% uzk(na,na) = real(uzk(na,na+1));
