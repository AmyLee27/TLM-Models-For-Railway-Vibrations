clear all
clc
close all

% addpath '.\asphalt_track_toolbox_v12'

warning('off','all')
warning

dofs = 3; % 1: 1 dof; 2: 2 dofs; 3: 3 dofs

if(1) %computation?
    %% PROPERTIES
    
    %TRACK DATA
    %--------------------
    width_t=1.5/2;                                %<-- track: gauge [m]
    
    bl_s=10;                                         %<-- substrate: width (bottom) [m]
    h_s=2;                                       %<-- substrate: height [m]
%     slope_s=30;
%     bu_s=bl_s-2*h_s/tand(slope_s);                  %<-- substrate: width (upper) [m]
    bu_s=4;
    z0_s=0;                                         %<-- substrate: z coordinate at bottom [m]
    n1_s=20;%12;                                        %<-- substrate: number horizontal elements (mesh) [pair]
    n2_s=20;                                         %<-- substrate: number vertical elements (mesh)
    Ezz_s=54E6;                                    %<-- substrate: Young's modulus [N/m^2]
    nu_s=0.35;                                       %<-- substrate: Poisson's ratio [-]
    rho_s=2000;                                     %<-- substrate: density [kg/m^3]
    D_s=0.03;                                          %<-- substrate: damping ratio
    Ezz_s=Ezz_s*(1+sqrt(-1)*D_s);                   %<-- substrate: complex stiffness [N/m^2]
    
    
    bl_a=3.5;                                         %<-- substrate: width (bottom) [m]
    bu_a=3.5;
    n1_a=round(n1_s*bl_a/bu_s);                    %<-- asphalt: number horizontal elements (mesh) [pair]
    if(mod(n1_a,2));n1_a=n1_a+1;end
    bl_a=bu_s/n1_s*n1_a;
    bu_a=bl_a;
    Ezz_a=2E9;                                      %<-- asphalt: Young's modulus [N/m^2]
    nu_a=0.35;                                      %<-- asphalt: Poisson's ratio [-]
    rho_a=2400;                                     %<-- asphalt: density [kg/m^3]
    D_a=0.015;                                       %<-- asphalt: damping ratio
    Ezz_a=Ezz_a*(1+sqrt(-1)*D_a);                   %<-- asphalt: complex stiffness [N/m^2]
    h_a=0.2;                                       %<-- asphalt: height [m]
    n2_a=2;                                         %<-- substrate: number vertical elements (mesh)
    
    h_b=0.20;                                       %<-- ballast: height [m]
    n2_b=2;                                         %<-- ballast: number vertical elements (mesh)
    % slope_b=30;
    % bl_b=bu_s;                                      %<-- ballast: width (bottom) [m]
    % bu_b=bl_b-2*h_b/tand(slope_b);                  %<-- ballast: width (upper) [m]
    bl_b=3;                                      %<-- ballast: width (bottom) [m]
    
    n1_b=round(n1_a*bl_b/bu_a);
    if(mod(n1_b,2));n1_b=n1_b+1;end
    bl_b=bu_a/n1_a*n1_b;
    bu_b=bl_b;
    
    Ezz_b=150E6;                                    %<-- ballast: Young's modulus [N/m^2]
    nu_b=0.25;                                       %<-- ballast: Poisson's ratio [-]
    rho_b=1800;                                     %<-- ballast: density [kg/m^3]
    D_b=0.05;                                          %<-- ballast: damping ratio
    Ezz_b=Ezz_b*(1+sqrt(-1)*D_b);                   %<-- ballast: complex stiffness [N/m^2]
    
    d_sl=0.65;                                         %<-- sleeper: distance [m]
    Ezz_sl=30E9;                                      %<-- sleeper: Young's modulus [N/m^2]
    b_sl=0.25;
    l_sl=2.4;
    D_sl=0.05;
    Ezz_sl=Ezz_sl*b_sl/d_sl*(1+sqrt(-1)*D_sl);
    nu_sl=0.2;                                        %<-- sleeper: Poisson's ratio [-]
    h_sl=0.2;                                       %<-- sleeper: height [m]
    rho_sl=2500;                                       %<-- sleeper: height [m]
    rho_sl=rho_sl*b_sl/d_sl;
    n2_sl=1;                                         %<-- sleeper: number vertical elements (mesh)
    
    n1_sl=round(n1_b*l_sl/bu_b);
    if(mod(n1_sl,2));n1_sl=n1_sl+1;end
    l_sl=bu_b/n1_b*n1_sl;
    
    Ezz_usl=30E9;                                     %<-- under sleeper pad: Young's modulus [N/m^2]
    nu_usl=0.2;                                        %<--  under sleeper pad: Poisson's ratio [-]
    rho_usl=2500;                                       %<--  under sleeper pad: density [kg/m^3]
    b_usl=0.25;
    l_usl=2.4;
    h_usl=0;                                       %<-- under sleeper pad: height [m]
    
    Ezz_usl=Ezz_usl*b_usl/d_sl;
    D_usl=0.05;                                       %<-- under sleeper pad: damping ratio
    Ezz_usl=Ezz_usl*(1+sqrt(-1)*D_usl);                   %<-- under sleeper pad: complex stiffness [N/m^2]
    
    rho_usl=rho_usl*b_usl/d_sl;
    n2_usl=1;                                         %<--  under sleeper pad: number vertical elements (mesh)
    
    n1_usl=round(n1_b*l_usl/bu_b);
    if(mod(n1_usl,2));n1_usl=n1_usl+1;end
    l_usl=bu_b/n1_b*n1_usl;
    
    kzz_rp=200E6;                                    %[N/m]
    b_rp=0.25;
    l_rp=0.1;
    h_rp=10e-3;                                       %<-- rail pad: height [m]
    
    Ezz_rp=80e6; %kzz_rp*h_rp/l_rp/b_rp;                                   %<-- rail pad: Young's modulus [N/m^2]
    nu_rp=0.48;                                        %<--  rail pad: Poisson's ratio [-]
    rho_rp=1300*b_rp/d_sl;                                       %<--  rail pad: density [kg/m^3]
    
    D_rp=0.08;                                    %<-- rail pad: damping [Ns^2/m]
    Ezz_rp=Ezz_rp*(1+sqrt(-1)*D_rp);                   %<-- rail pad: complex stiffness [N/m^2]
    
    Ezz_rp=Ezz_rp*b_rp/d_sl;                                   %<-- rail pad: stiffness [N/m/m]
    kzz_rp=kzz_rp*b_rp/d_sl;                                   %<-- rail pad: stiffness [N/m/m]
    
    h_r=0.155;
    l_r=l_rp;                
    A_r=0.0155;                                    %<-- rail: Section [m^2]
    Ixx_r=3.103e-5;                                %<-- rail: bending stiffness [m^4]-Ix
    Izz_r=1.292e-5;                                 % Iy
    E_r=2.1e11;                                     %<-- rail: Young's modulus [N/m^2]
    nu_r=0.1;                                       %<--  rail: Poisson's ratio [-]
    rho_r=7850;                                     %<-- rail: Density [kg/m^3]
    D_r=0.08;                                       %<-- rail: damping ratio
    E_r=E_r*(1+sqrt(-1)*D_r);                       %<-- rail: complex stiffness [N/m^2]
    rail_type = 1;                                  % rail type: 1-solid element; 2-beam element
    
    
    % SOIL DATA
    %------------------------
    % Soil density
    var1_d = {'p'};
    var1_v = [2000]';
    % Secondary wave speed
    var2_d = {'Vs'};
    var2_v = [100];
    % Primary wave speed
    var3_d = {'v'};
    var3_v = [0.35]';
    % Soil layer heights
    layer_heights =[ 0]';
    % Dampings of soil in each layer
    damping = [0.03 ]';
    
    for x=1:length(var1_v)
        [E(x) G(x) K(x) lambda(x) M(x) Poisson(x) Density(x) Vs(x) Vp(x) Vr(x)] = elastodynamic(var1_d,var2_d,var3_d,var1_v(x),var2_v(x),var3_v(x));
    end
    
    i = sqrt(-1);
    
    [soil_prop_upper,soil_prop_inf,element_number,elements,element_coor,cinf] = soil_prop_multi(layer_heights, damping,var1_d,var2_d,var3_d,var1_v,var2_v,var3_v);
    soil_prop_upper = [soil_prop_upper(:,1) flipud(soil_prop_upper(:,2:5)) soil_prop_upper(:,1)];
    nodes1 = [elements element_coor];  % nodes information
    elements1 = element_number;
    
    [feNod,feElt,feTyp,feMat,feSec,P]=asphalt_track_8nodes_test(width_t,bl_s,h_s,bu_s,z0_s,n1_s,n2_s,h_b,bu_b,bl_b,n2_b,n1_b,n2_a,n1_a,h_a,rho_s,Ezz_s,nu_s,...
        rho_b,Ezz_b,nu_b,rho_a,Ezz_a,nu_a,bl_a,bu_a,rho_usl,Ezz_usl,nu_usl,h_usl,n2_usl,n1_usl,l_usl,rho_sl,Ezz_sl,nu_sl,h_sl,n2_sl,n1_sl,l_sl,...
        Ezz_rp,rho_rp,nu_rp,h_rp,l_rp,E_r,l_r,h_r,rho_r,nu_r,A_r,Ixx_r,Izz_r,dofs,rail_type);
    
    [M,K0,K1,K2,K3,K4,C0,feDof]=k(feNod,feElt,feTyp,feMat,feSec,dofs);
    nFeDof=size(M,1);
    nFeElt=size(feElt,1);
    
    figure
    hold on
    X=zeros(8,nFeElt);
    Y=zeros(8,nFeElt);
    for ielem=1:nFeElt
        if(length(nonzeros(feElt(ielem,5:12)))==8)
            X(:,ielem)=feNod(feElt(ielem,5:12),2);
            Y(:,ielem)=feNod(feElt(ielem,5:12),3);
            plot(X(1:8,ielem),Y(1:8,ielem),'ob')
        else
            X(1:7,ielem)=feNod(feElt(ielem,5:11),2);
            Y(1:7,ielem)=feNod(feElt(ielem,5:11),3);
            plot(X(1:7,ielem),Y(1:7,ielem),'xr')
        end
    
    end
    
    k1=0;
    k2=-20:10/1024:20;
    c=0;
    nWave=length(k1);
    f = 20;
    
    zeroNod=find(feNod(:,3)==0); % find the index of the bottom nodes
    nEle = (length(zeroNod)-1)/2;
    nitmnode = zeros(nEle,1);
    % find the coordinates of the middle nodes at the bottom
    for k=1:nEle
        nitmnode(k,1) = feNod(zeroNod(2*k),2) ;
    end
    
    
    % Calculate green's function in each direction
%     [uxxk,uyxk,uzxk]=greenfunction_PAC_test(k1,k2,cinf,f,c,nodes1,elements1,soil_prop_upper,soil_prop_inf,1);
%     [uxyk,uyyk,uzyk]=greenfunction_PAC_test(k1,k2,cinf,f,c,nodes1,elements1,soil_prop_upper,soil_prop_inf,2);
%     [uxzk,uyzk,uzzk]=greenfunction_PAC_test(k1,k2,cinf,f,c,nodes1,elements1,soil_prop_upper,soil_prop_inf,3);
% 
%     figure
%     plot(k2,uxxk, k2, uyyk, k2, uzyk, k2, uzzk)
%     legend('xx','yy','zy','zz')
    
%     figure
%     plot(k2,uxxk2, k2, uyyk2, k2, uzyk2, k2, uzzk2)
%     legend('xx','yy','zy','zz')
% %     
    zeroNod=find(feNod(:,3)==0);
    u=zeros(nFeDof,nWave);
    duy=zeros(nFeDof,nWave);
    u2=zeros(nFeDof,1);
    
    for iWave=1:nWave
            ky=k1(iWave);
            jk1=iWave;
            omega = 2*pi*f -ky*c;
            
%             [Kdyn,G]=KdynMatrix_PAC(uxxk,uyxk,uzxk,uxyk,uyyk,uzyk,uxzk,uyzk,uzzk,k2,jk1,nitmnode,zeroNod,feNod);
            %FEM Matrices
            Kt=K0 + sqrt(-1)*omega*C0 - sqrt(-1)*K1*ky - K2*ky^2 + sqrt(-1)*K3*ky^3 + K4*ky^4 - omega^2*M;
%             Kt=T1*Kt+t1;

            
            %Soil stiffness
            Kts=zeros(size(Kt,1),size(Kt,2));
            Kts(1:size(Kt,1),1:size(Kt,2))=Kt;
            
            % add the stiffness matrix to the global matrix
%             KtEq = zeros(size(Kt,1),size(Kt,2));
%             KtEq(1:numel(zeroNod)*dofs,1:numel(zeroNod)*dofs) = Kdyn;
%             Kts = Kts+KtEq;

            % for the case it is rigid bottom
            Kts = Kdyn_update(Kts,zeroNod);

            Kts = sparse(Kts);
            
            u2=Kts\P;
%             u(1:1:end-3,iWave,iFreq)=u2;% depends on dofs
            u(1:1:end,iWave) = u2;
%             duy(1:1:end-3,iWave,iFreq)=-sqrt(-1)*ky*u2;
            duy(1:1:end,iWave)=-sqrt(-1)*ky*u2;
    end
      
    figure
    X = feNod(:,2);
    Y = feNod(:,3);
    tri=delaunay(X,Y);
    trisurf(tri,X,Y,abs(squeeze(u(3:3:end,1))));
    shading('interp')
    view([0 -90])
    hold on
    plot(X,Y,'ok')
    title('Embankment - 5Hz - KD')
    set(gca,'Ydir','reverse')
    hold off
    c=colorbar;
    drawnow    
    
    figure
    X = nodes(:,2);
    Y = nodes(:,3);
    tri=delaunay(X,Y);
    trisurf(tri,X,Y,abs(squeeze(uz4)));
    shading('interp')
    view([0 -90])
    hold on
    plot(X,Y,'ok')
    title('Embankment - 5Hz - PAC')
    set(gca,'Ydir','reverse')
    hold off
    c=colorbar;
    drawnow 
    
    figure
    plot(1:1281, abs(u(3:3:end,1)), 1:1281, abs(uz4))
    legend('KD', 'PAC')
    title('Embankment - 5Hz - all nodes comparison')
    
      figure
    plot(1:1281, u(3:3:end,1), 1:1281, uz4)
    legend('KD', 'PAC')
    title('Slab - 5Hz - all nodes comparison')
    
end  
    
    
    
    
    
    
    
    
    
    
    