function [uzk,strt1,strt2,strt3,strt4,strt5,strt6]=greenfunction_strain(k1max,k1inc,cinf,v,nodes1,elements1,properties1,properinf1)

ki=-k1max:k1inc:k1max;   % define the wavenumbers in -x and -y directions
kj=ki;

elestrain = 1:length(elements1(:,1));
uzk=zeros(length(ki),length(kj)); % disp infor for the top node

% Set up matrices for strain values (6 in total, xx,yy,zz,xy,yz,xz)
strt1=zeros(length(ki),length(ki),length(elestrain));
strt2=zeros(length(ki),length(ki),length(elestrain));
strt3=zeros(length(ki),length(ki),length(elestrain));
strt4=zeros(length(ki),length(ki),length(elestrain));
strt5=zeros(length(ki),length(ki),length(elestrain));
strt6=zeros(length(ki),length(ki),length(elestrain));

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
    [uzk2,str1,str2,str3,str4,str5,str6]=Mnoddispdiv1_strain(kj,k1,w,nodes1,cinf,elements1,properinf1,nodinf1,K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F,no,elestrain,properties1);
    uzk(j,:)=uzk2;
    
    % Store the strain information in all directions (k1.k2,z)
    strt1(j,:,:)=str1;
    strt2(j,:,:)=str2;
    strt3(j,:,:)=str3;
    strt4(j,:,:)=str4;
    strt5(j,:,:)=str5;
    strt6(j,:,:)=str6;
       
end

% Get rid of NaN values in the disp, strain, stress results
na = find(ki==0);
uzk(na,na) = real(uzk(na,na+1));

for i = 1:length(elestrain)
    strt1(na,na,i) = real(strt1(na,na+1,i));
    strt2(na,na,i) = real(strt2(na,na+1,i));
    strt3(na,na,i) = real(strt3(na,na+1,i));
    strt4(na,na,i) = real(strt4(na,na+1,i));
    strt5(na,na,i) = real(strt5(na,na+1,i));
    strt6(na,na,i) = real(strt6(na,na+1,i));
end

