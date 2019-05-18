%%%%%%Subrotine to compute the strain and the stress%%%%%%%%
%%%%%Developed by Pedro Alves Costa - 06/12/2007-30/09/2008

function [strt1,strt2,strt3,strt4,strt5,strt6]=Strain_p(ax,ay,az,elements,nodes,properties,k1,k2,elestrain)

%%%elestrain - number of the elements where the strains and stresses are
%%%computed
%%%code - 1 - computes the strains on the gauss points
% 2 - computes the strains on the center of the element

Gauss1;
PG1(:,1)=[];


 for j=1:length(elestrain(1,:)),
     
    axi=ax(elements(elestrain(j),:),1);
    ayi=ay(elements(elestrain(j),:),1);
    azi=az(elements(elestrain(j),:),1);
     
    for s=1:3, %Organizes the vector of displacemnets of the nodes of the element j
         u(3*s-2,1)=axi(s,1);
         u(3*s-1,1)=ayi(s,1);
         u(3*s,1)=azi(s,1);
    end
%     
%     if code==1,
%         for s=1:size(PG1,1),
%             Y=nodes(elements(elestrain(j),:),:);  %coordinates of the nodes of element j
%             [B]=Bmatrix(Y,PG1(s,1),PG1(s,2),k);
%             stre(s,:)=(B*u);               %Computes the strain
%             elety=properties(elestrain(j),1);
%             E=properties(elestrain(j),2);
%             niu=properties(elestrain(j),3);
%             del=properties(elestrain(j),4);
%             Et=properties(elestrain(j),6);
%             niut=properties(elestrain(j),7);
%             Gt=properties(elestrain(j),8);
%             stres(s,:)=(Dmatrix(E,niu,del,w,Et,niut,Gt,elety)*(stre(s,:)).').';    %computes the stress
%         end
%         strt(4*j-3:4*j,:)=stre;
%         strest(4*j-3:4*j,:)=stres;
%     else
         Y=nodes(elements(elestrain(j),:),:);  %coordinates of the nodes of element j
         [B]=Bmatrix1_p(Y,0,k1,k2);
         stre(1,:)=(B*u);               %Computes the strain
 %       elety=properties(elestrain(j),1);
% %         E=properties(elestrain(j),2);
% %         niu=properties(elestrain(j),3);
% %         del=properties(elestrain(j),4);
% %         Et=properties(elestrain(j),6);
% %         niut=properties(elestrain(j),7);
% %         Gt=properties(elestrain(j),8);
% %         stres(1,:)=(Dmatrix(E,niu,del,w,Et,niut,Gt,elety)*(stre(1,:)).').';    %computes the stress
% %      

% It is worth nothing that the strains obtained now are in an inverse order,
% i.e., from bottom to top

        strt1(j)=stre(1,1);
        strt2(j)=stre(1,2);
        strt3(j)=stre(1,3);
        strt4(j)=stre(1,4);
        strt5(j)=stre(1,5);
        strt6(j)=stre(1,6);
        
 end

