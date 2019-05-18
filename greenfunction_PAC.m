function [uxk,uyk,uzk]=greenfunction_PAC(k1max,k1inc,cinf,v,nodes1,elements1,properties1,properinf1,dir)

ki=-k1max:k1inc:k1max;   % define the wavenumbers in -x and -y directions
ki(ki==0) = 0.01;
kj=ki;

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
    w=-k1*v; % k(wavenumber) = 2*pi*f/v = w(angular frequency)/v;
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
% na = find(ki==0);
% uzk(na,na) = real(uzk(na,na+1));
