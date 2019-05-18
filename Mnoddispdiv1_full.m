function [uzk2,str1,str2,str3,str4,str5,str6,stres1,stres2,stres3,stres4,stres5,stres6]=Mnoddispdiv1_full(kj,k1,w,nodes1,cinf,elements1,properinf1,nodinf1,K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F,no,elestrain,properties1)

   for s=1:length(kj);
        k2=kj(s);
        [ax,ay,az]=noddispdiv1(k1,k2,w,nodes1,cinf,properinf1,nodinf1,K0f,K1f,K2f,K3f,K4f,K5f,Mf,Kr,F); % obtain disp in all depths 
        [uzk2(s,:),uyk2(s,:),uxk2(s,:)]=forforfor(no,az,ay,ax); % obtain the disp on the top layer
        [strt1,strt2,strt3,strt4,strt5,strt6]=Strain_p(ax,ay,az,elements1,nodes1,properinf1,k1,k2,elestrain);
        [stre1,stre2,stre3,stre4,stre5,stre6]=Stress_p(strt1,strt2,strt3,strt4,strt5,strt6,properties1,w,elestrain);
        
        % rearrage the strains from top to bottom, and store them in rows
        % i.e., in k2 direction, which is transversal direction (-y)
        str1(s,:) = strt1;   % size of strain or stress is (k1, length(elastrain)), i.e., k by z
        str2(s,:) = strt2;
        str3(s,:) = strt3;
        str4(s,:) = strt4;
        str5(s,:) = strt5;
        str6(s,:) = strt6;
        
        stres1(s,:) = stre1;   
        stres2(s,:) = stre2;
        stres3(s,:) = stre3;
        stres4(s,:) = stre4;
        stres5(s,:) = stre5;
        stres6(s,:) = stre6;

   end