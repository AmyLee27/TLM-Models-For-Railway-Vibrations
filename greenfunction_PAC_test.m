function [uxk,uyk,uzk]=greenfunction_PAC_test(k1,k2,cinf,f,v,nodes1,elements1,properties1,properinf1,dir)

ki = k1;
kj=k2;

elestrain = 1:length(elements1(:,1));
uzk=zeros(length(ki),length(kj)); % disp infor for the top node
uyk=zeros(length(ki),length(kj));
uxk=zeros(length(ki),length(kj));


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

switch dir
    case 1  % force on x-direction
        loads1=[length(nodes1(:,1)),1,0,0];
    case 2  % force on y-direction
        loads1=[length(nodes1(:,1)),0,1,0];
    case 3  % force on z-direction
        loads1=[length(nodes1(:,1)), 0,0,1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F]=MainMatrix1(nodes1,elements1,properties1,suports1,loads1);

parfor j=1:length(ki);
    
    k1=ki(j);  % It is calculated in k1 direction, i.e., longitudinal direction (-x)
    w=2*pi*f - k1*v; % k(wavenumber) = 2*pi*f/v = w(angular frequency)/v;
   [uzk2,uyk2,uxk2]=Mnoddispdiv1_PAC(kj,k1,w,nodes1,cinf,elements1,properinf1,nodinf1,K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F,no,elestrain,properties1);

    uzk(j,:)=uzk2;
    uyk(j,:)=uyk2;
    uxk(j,:)=uxk2;

end
% uzk(length(ki),length(kj))=uzk(length(ki),length(kj)-1);
% [US]=Reordek2(uzk,kj);
% [US]=Reorde(US,ki);
% uzk=US;
% 
na = find(kj==0);
uzk(1,na) = uzk(1,na+1);
uyk(1,na) = uyk(1,na+1);
uxk(1,na) = uxk(1,na+1);
