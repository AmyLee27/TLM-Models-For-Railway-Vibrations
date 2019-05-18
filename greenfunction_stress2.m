function [uzk,stress1,stress2,stress3,stress4,stress5,stress6]=greenfunction_stress2(knum,cinf,v,nodes1,elements1,properties1,properinf1)

ki=-knum;
kj=ki;

elestrain = 1:length(elements1(:,1));
uzk=zeros(length(ki),length(kj)); % disp infor for the top node

% Set up matrices for stress values (6 in total, xx,yy,zz,xy,yz,xz)
stress1=zeros(length(ki),length(ki),length(elestrain));
stress2=zeros(length(ki),length(ki),length(elestrain));
stress3=zeros(length(ki),length(ki),length(elestrain));
stress4=zeros(length(ki),length(ki),length(elestrain));
stress5=zeros(length(ki),length(ki),length(elestrain));
stress6=zeros(length(ki),length(ki),length(elestrain));

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
%     [uzk2,dispx,dispy,dispz,str1,str2,str3,str4,str5,str6,stres1,stres2,stres3,stres4,stres5,stres6]=Mnoddispdiv1(kj,k1,w,nodes1,cinf,elements1,properinf1,nodinf1,K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F,no,elestrain,properties1);
    [uzk2,stres1,stres2,stres3,stres4,stres5,stres6]=Mnoddispdiv1_stress(kj,k1,w,nodes1,cinf,elements1,properinf1,nodinf1,K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F,no,elestrain,properties1);
    
    uzk(j,:)=uzk2;    

    % Store the stress information in all directions (k1.k2,z)
    stress1(j,:,:)=stres1;
    stress2(j,:,:)=stres2;
    stress3(j,:,:)=stres3;
    stress4(j,:,:)=stres4;
    stress5(j,:,:)=stres5;
    stress6(j,:,:)=stres6;
           
end

% Get rid of NaN values in the disp, strain, stress results
na = find(ki==0);
uzk(na,na) = real(uzk(na,na+1));

for i = 1:length(elestrain)
    stress1(na,na,i) = real(stress1(na,na+1,i));
    stress2(na,na,i) = real(stress2(na,na+1,i));
    stress3(na,na,i) = real(stress3(na,na+1,i));
    stress4(na,na,i) = real(stress4(na,na+1,i));
    stress5(na,na,i) = real(stress5(na,na+1,i));
    stress6(na,na,i) = real(stress6(na,na+1,i));
end



