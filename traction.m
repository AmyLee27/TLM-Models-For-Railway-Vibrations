function [traction,uf,duf] = traction(nWave,nFreq,nFeDof,k1,k2,f,c,uxxk,uyxk,uzxk,uxyk,uyyk,uzyk,uxzk,uyzk,uzzk,nitmnode,zeroNod,feNod,midDof,dofs)

u=zeros(nFeDof,nWave);
duy=zeros(nFeDof,nWave);

unew = zeros(nFeDof,nWave,nFreq);
duynew=zeros(nFeDof,nWave,nFreq);

uf=zeros(nFeDof,nFreq);

for iFreq=1:nFreq

for iWave=1:nWave
        
        ky=k1(iWave);
        jk1=iWave;
        omega = 2*pi*f -ky*c;
            
        [Kdyn,G]=KdynMatrix_PAC(uxxk,uyxk,uzxk,uxyk,uyyk,uzyk,uxzk,uyzk,uzzk,k2,jk1,nitmnode,zeroNod,feNod);
        % FEM Matrices
        Kt=K0 + sqrt(-1)*omega*C0 - sqrt(-1)*K1*ky - K2*ky^2 + sqrt(-1)*K3*ky^3 + K4*ky^4 - omega^2*M;
%       Kt=T1*Kt+t1;

            
        % Soil stiffness
        Kts=zeros(size(Kt,1),size(Kt,2));
        Kts(1:size(Kt,1),1:size(Kt,2))=Kt;
            
        % add the stiffness matrix to the global matrix
        KtEq = zeros(size(Kt,1),size(Kt,2));
        KtEq(1:numel(zeroNod)*dofs,1:numel(zeroNod)*dofs) = Kdyn;
        Kts = Kts+KtEq;

        Kts = sparse(Kts);
            
        u2=Kts\P;
%       u(1:1:end-3,iWave,iFreq)=u2;% depends on dofs
        u(1:1:end,iWave) = u2; % node disp
%       duy(1:1:end-3,iWave,iFreq)=-sqrt(-1)*ky*u2;
        duy(1:1:end,iWave) = -sqrt(-1)*ky*u2;
        
        % calculate the traction at the middle node of the bottom track for all DOFs       
        traction(:,iWave) = G\u(midDof,iWave);       

        % calculate free field response
%         uf(iWave) = field_PAC(uzzk,k2,jk1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,iWave));

end
unew(:,:,iFreq) = u;% dispalcement in all dofs (include rail, railpad, etc etc) all freqs all wavenumbers
duynew(:,:,iFreq) = duy; % duy = du/dy - derivative
end

%FREQUENCY RESPONSE
    y=0.001;
    for iDof=1:nFeDof
        for iFreq=1:nFreq
            ky=k1;
            uf(iDof,iFreq)=1/pi*intfilon(ky,-sqrt(-1)*y,unew(iDof,:,iFreq)); % transform to frequency domain only
            duf(iDof,iFreq)=1/pi*intfilon(ky,-sqrt(-1)*y,duynew(iDof,:,iFreq));
        end
    end