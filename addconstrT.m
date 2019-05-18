function [T,t]=addconstrT(Constr,T,t)

dof=Constr(:,3:2:size(Constr,2));
Coef=-Constr(:,2:2:(size(Constr,2)-1))./repmat(Constr(:,2),1,size(dof,2));

SlaveDOF=dof(:,1);
nEqn=size(Constr,1);

nDOF=size(T,1);

T=eye(nDOF);
t=zeros(nDOF);
DOF=1:1:nDOF;

for iEqn=1:nEqn
    locSlave=find(DOF==SlaveDOF(iEqn,1));
    if size(dof,2)>1
        for idof=2:size(dof,2)
            locdof=find(DOF==dof(iEqn,idof));
            T(locdof,locSlave)=Coef(iEqn,idof);
            T(locdof,locdof)=1;
            T(locSlave,locSlave)=0;
        end
    end
    if size(dof,2)>1
        for idof=1:size(dof,2)
            locdof=find(DOF==dof(iEqn,idof));
            t(locSlave,locdof)=-Coef(iEqn,idof);
        end
    end
end
