function [uzk]=vibfemgreendiv(k1max,k1inc,cinf,v)
%cinf - código - infinitos =1, caso contrário 0
%warning off;
ki=-k1max:k1inc:0;
kj=ki;
%kj=-k2max:k2inc:0;

uzk=zeros(length(ki),length(kj));


%Finite elements properties
 load nodes1.dat;
 no=max(nodes1(:,1));
 nodes1(:,1)=[];
 load elements1.dat;
 elements1(:,1)=[];
 load properties1.dat;
 properties1(:,1)=[];
 %%%%%%%%%%%%%%%%%%%%%%%%%
 
 %Infinite elements properties
 if cinf==1;
     load properinf1.dat;
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
    k1=ki(j);
    w=-k1*v;
    [uzk2]=Mnoddispdiv1(kj,k1,w,nodes1,cinf,properinf1,nodinf1,K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F,no);  
    uzk(j,:)=uzk2;
    
end
uzk(length(ki),length(kj))=uzk(length(ki),length(kj)-1);
[US]=Reordek2(uzk,kj);
[US]=Reorde(US,ki);
uzk=US;


%uxk=Reorde2(uxk,ki);
%uyk=Reorde(uyk,ki);
%uzk=Reorde(uzk,ki);
