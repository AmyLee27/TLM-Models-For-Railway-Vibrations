function dof=femdof(elt,typ,dofs)

nElt=size(elt,1);
for iElt=1:nElt
  eltdef2=eltdef(typ,elt(iElt,2),dofs);
  dof=zeros([eltdef2.nDof,1]);
  for iDof=1:eltdef2.nDof
    eltnod=floor(eltdef2.Dof(iDof));
    dof(iDof)=elt(iElt,4+eltnod)+(eltdef2.Dof(iDof)-eltnod);
  end
  eltdof{iElt}=dof(:).';
end
dof=sort(unique(round(cell2mat(eltdof)*100)))/100;
