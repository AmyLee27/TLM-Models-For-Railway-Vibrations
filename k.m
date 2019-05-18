function [M,K0,K1,K2,K3,K4,C0,dof]=k(nod,elt,typ,mat,sec,dofs)

dof=femdof(elt,typ,dofs); 

nDof=numel(dof);
nElt=size(elt,1);

for iElt=1:nElt
  def=eltdef(typ,elt(iElt,2),dofs);
  eltdofind=zeros([1,def.nDof]);
  for iDof=1:def.nDof
    eltnod=floor(def.Dof(iDof));
    idof=def.Dof(iDof)-eltnod;
    util=find(abs(dof-(elt(iElt,4+eltnod)+idof))<0.0001,1); 
    if isempty(util)
      eltdofind(iDof)=0;
    else
      eltdofind(iDof)=util;
    end
  end
  utilind=find(eltdofind);
  nutilind=length(utilind);
  iind{iElt}=repmat(eltdofind(utilind),1,nutilind) ;
  jind{iElt}=reshape(repmat(eltdofind(utilind),nutilind,1),1,nutilind^2);
  
  [Me,K0e,K1e,K2e,K3e,K4e,C0e]=ke(nod,elt(iElt,:),typ,mat(iElt,2:5),sec,dofs);
 
  sM{iElt} =reshape(Me(utilind,utilind),1,nutilind^2);
  sK0{iElt}=reshape(K0e(utilind,utilind),1,nutilind^2);
  sK1{iElt}=reshape(K1e(utilind,utilind),1,nutilind^2);
  sK2{iElt}=reshape(K2e(utilind,utilind),1,nutilind^2);
  sK3{iElt}=reshape(K3e(utilind,utilind),1,nutilind^2);
  sK4{iElt}=reshape(K4e(utilind,utilind),1,nutilind^2);
  sC0{iElt}=reshape(C0e(utilind,utilind),1,nutilind^2);
end


iind=cell2mat(iind);
jind=cell2mat(jind);
ind=sub2ind([nDof nDof],iind,jind);
[dum,colsort]=sort(ind);

sM=cell2mat(sM);
sK0=cell2mat(sK0);
sK1=cell2mat(sK1);
sK2=cell2mat(sK2);
sK3=cell2mat(sK3);
sK4=cell2mat(sK4);
sC0=cell2mat(sC0);

M = sparse(iind(colsort),jind(colsort),sM(colsort),nDof,nDof);
K0= sparse(iind(colsort),jind(colsort),sK0(colsort),nDof,nDof);
K1= sparse(iind(colsort),jind(colsort),sK1(colsort),nDof,nDof);
K2= sparse(iind(colsort),jind(colsort),sK2(colsort),nDof,nDof);
K3= sparse(iind(colsort),jind(colsort),sK3(colsort),nDof,nDof);
K4= sparse(iind(colsort),jind(colsort),sK4(colsort),nDof,nDof);
C0= sparse(iind(colsort),jind(colsort),sC0(colsort),nDof,nDof);

