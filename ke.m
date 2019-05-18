function [Me,K0e,K1e,K2e,K3e,K4e,C0e]=ke(nod,elt,typ,mat,sec,dofs)

nElt=size(elt,1);
nNod=size(nod,1);
nTyp=size(typ,1);

typID=elt(2);
matID=elt(3);
secID=elt(4);
eltdef2=eltdef(typ,typID,dofs);

nodcoord=zeros(eltdef2.nNod,2);
for iNod=1:eltdef2.nNod
    nodcoord(iNod,:)=id2prop(nod,elt(4+iNod));
end

typName = elttypname(typ,typID);

switch typName
    case 'volu4'
        matutil=id2prop(mat,matID);
        E=matutil(1);
        nu=matutil(2);
        rho=matutil(3);
        [Me,K0e,K1e,K2e,K3e,K4e,C0e]=volu4(nodcoord,E,nu,rho,dofs);
    case 'volu7'
        matutil=id2prop(mat,matID);
        E=matutil(1);
        nu=matutil(2);
        rho=matutil(3);
        [Me,K0e,K1e,K2e,K3e,K4e,C0e]=volu7(nodcoord,E,nu,rho,dofs);
    case 'volu8'
        matutil=id2prop(mat,matID);
        E=matutil(1);
        nu=matutil(2);
        rho=matutil(3);
%         [Me,K0e,K1e,K2e,K3e,K4e,C0e]=volu8(nodcoord,E,nu,rho,dofs); 
        [Me,K0e,K1e,K2e,K3e,K4e,C0e]=volu8_PAC(nodcoord,E,nu,rho,dofs); 
    case 'beam1'
        matutil=id2prop(mat,matID);
        E=matutil(1);
        nu=matutil(2);
        rho=matutil(3);
        secutil=id2prop(sec,secID);
        A=secutil(1);
        Ixx=secutil(2);
        Izz=secutil(3);
        [Me,K0e,K1e,K2e,K3e,K4e,C0e]=beam1(E,nu,rho,A,Ixx,Izz,dofs);
    case 'mass1'
        secutil=id2prop(sec,secID);
        m=secutil(1);
        [Me,K0e,K1e,K2e,K3e,K4e,C0e]=mass1(m);
    case 'link2'  % 2 node 2.5D link
        secutil=id2prop(sec,secID);
        kzz=secutil(1);
        czz=secutil(2);
        [Me,K0e,K1e,K2e,K3e,K4e,C0e]=link2(kzz,czz);
end

%===============================================================================
function [prop,propind]=id2prop(item,itemID)

nID=numel(itemID);
propind=[];

for iID=1:nID
    ind=find(item(:,1)==itemID(iID),1);
    propind=[propind,ind];
end
prop=item(propind,2:end);

%===============================================================================
function typName = elttypname(typ,typID)
nTyp=size(typ,1);
typName='';
for iTyp=1:nTyp
    if typ{iTyp,1}==typID, typName=typ{iTyp,2}; end
end

%===============================================================================
function opt = typ2opt(typ,typID)
opt={};
nTyp=size(typ,1);
for iTyp=1:nTyp
    if ((typ{iTyp,1}==typID) && (length(typ(1,:))==3))
        opt=typ{iTyp,3};
    end
end
