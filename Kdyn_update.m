function Kdyn = Kdyn_update(Kdyn, zeroNode)

% replace the stiffness matrix of FEM to set boundary condition

nLen = length(zeroNode);
for j=1:nLen
    
    Kdyn((3*j-2):(3*j),:) = 0;
    Kdyn(3*j-2, 3*j-2) = 1;
    Kdyn(3*j-1, 3*j-1) = 1;
    Kdyn(3*j, 3*j) = 1;

end