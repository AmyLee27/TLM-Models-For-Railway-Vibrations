function [epsilonn, sigman] = track_element(nFeDof,nFeElt,dofs,nFreq,Ezz_s,nu_s,Ezz_a,nu_a,Ezz_b,nu_b,Ezz_usl,nu_usl,Ezz_sl,nu_sl,Ezz_rp,nu_rp,E_r,nu_r,uf,duf)

%STRAIN STRESS SOLID ELEMENT RESULTS
    
    % 6 strain/stress elements: xx,yy,zz,xy,zy,zx
    
    stress=1;
    sigman=zeros(6,nFeElt,nFreq); % final output
    sigman2=zeros(6,nFeElt,nFeDof/dofs,nFreq);
    sigma=zeros(6,8,nFeElt,nFreq);
    
    epsilonn=zeros(6,nFeElt,nFreq); % final output
    epsilonn2=zeros(6,nFeElt,nFeDof/dofs,nFreq);
    epsilon=zeros(6,8,nFeElt,nFreq);% first value 6 means the number of strain elements                                   
     
    
    if(stress)

        for iFreq=1:nFreq
            for iElt=1:nFeElt
                if(min(feElt(iElt,5:12))>0)
                    nNode=8;
                elseif(min(feElt(iElt,5:12))==0)
                    nNode=7;
                end
                nodesaux=feElt(iElt,5:(4+nNode));
                nodes=zeros(1,nNode*dofs);
                for inodes=1:nNode
                    nodes(1,(inodes-1)*dofs+[1:dofs])=(nodesaux(inodes)-1)*dofs+[1:dofs];%vertical nodes
                end
                
                mat=feElt(iElt,3);
                nodcoord=feNod(nodesaux,2:3);
                switch(mat)
                    case 1
                        E=Ezz_s;
                        nu=nu_s;
                    case 2
                        E=Ezz_a;
                        nu=nu_a;
                    case 3
                        E=Ezz_b;
                        nu=nu_b;
                    case 4
                        E=Ezz_usl;
                        nu=nu_usl;
                    case 5
                        E=Ezz_sl;
                        nu=nu_sl;
                    case 6
                        E=Ezz_rp;
                        nu=nu_rp;
                    case 7
                        E=E_r;
                        nu=nu_r;
                end
                if nNode==8
                    [epsilon(:,:,iElt,iFreq),sigma(:,:,iElt,iFreq)]=volu8s(nodcoord,uf(nodes,iFreq),duf(nodes,iFreq),E,nu);
                elseif nNode==7
                    [epsilon(:,1:nNode,iElt,iFreq),sigma(:,1:nNode,iElt,iFreq)]=volu7s(nodcoord,uf(nodes,iFreq),duf(nodes,iFreq),E,nu);
                end
            end
        end
        %write element results
        for iFreq=1:nFreq
            for iElt=1:1:nFeElt
                if(min(feElt(iElt,5:12))>0)
                    nodes=feElt(iElt,5:12);
                    epsilonn2(:,iElt,nodes,iFreq)=squeeze(epsilon(:,:,iElt,iFreq));
                    sigman2(:,iElt,nodes,iFreq)=squeeze(sigma(:,:,iElt,iFreq));
                    epsilonn(:,iElt,iFreq) = mean(epsilonn2(:,iElt,nodes,iFreq),2);
                    sigman(:,iElt,iFreq) = mean(sigman2(:,iElt,nodes,iFreq),2);
                elseif(min(feElt(iElt,5:12))==0)
                    nodes=feElt(iElt,5:11);
                    epsilonn2(:,iElt,nodes,iFreq)=squeeze(epsilon(:,1:7,iElt,iFreq));
                    sigman2(:,iElt,nodes,iFreq)=squeeze(sigma(:,1:7,iElt,iFreq));
                    epsilonn(:,iElt,iFreq) = mean(epsilonn2(:,iElt,nodes,iFreq),2);
                    sigman(:,iElt,iFreq) = mean(sigman2(:,iElt,nodes,iFreq),2);
                end
            end
        end
        
        
    end