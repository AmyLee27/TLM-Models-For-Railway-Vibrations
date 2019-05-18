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
    
    bl_s=4;                                         %<-- substrate: width (bottom) [m]
    h_s=0.3;                                       %<-- substrate: height [m]
%     slope_s=30;
%     bu_s=bl_s-2*h_s/tand(slope_s);                  %<-- substrate: width (upper) [m]
    bu_s=4;
    z0_s=0;                                         %<-- substrate: z coordinate at bottom [m]
    n1_s=20;%12;                                        %<-- substrate: number horizontal elements (mesh) [pair]
    n2_s=2;                                         %<-- substrate: number vertical elements (mesh)
    Ezz_s=30000E6;                                    %<-- substrate: Young's modulus [N/m^2]
    nu_s=0.35;                                       %<-- substrate: Poisson's ratio [-]
    rho_s=1900;                                     %<-- substrate: density [kg/m^3]
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
    var1_v = [1800]';
    % Secondary wave speed
    var2_d = {'G'};
    var2_v = [10E6];
    % Primary wave speed
    var3_d = {'v'};
    var3_v = [0.45]';
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
    soil_prop_upper_original = soil_prop_upper;
    soil_prop_inf_original = soil_prop_inf;
    nodes1 = [elements element_coor];  % nodes information
    elements1 = element_number;
    
    %--------------------------------------------------------------------------
    %   Define initial variables for non-linearity soil
    %--------------------------------------------------------------------------

    itmax = 7;    % maximum number of iterations
    tol = 0.03;    % convergence tolerance is 3%
    Ip = 0;  % Plasticity index - could be adjusted based on needs
    alpha = 0.65;   % Coefficient for octhedral strain

    track_weight = ( rho_s*h_s*bl_s )*9.81; % track weight including rail, sleeper and ballast

    % Input the matrices needed in the Green's function computation
    strain_depth = nodes1(elements1(:,4),2);  % layer depth for confining stress calculation (pick up the bottom node of each thin layer)
    sigma1(1) = (soil_prop_upper(1,5)*strain_depth(1)*9.81)/1000;
    for n = 2:length(soil_prop_upper(:,5))
        sigma1(n) = sigma1(n-1) + (soil_prop_upper(n,5)*(strain_depth(n)-strain_depth(n-1))*9.81)/1000;  % Confining stresses of each layer in KPa
    end
    sigma1 = fliplr(sigma1);
    poisson_ratio = soil_prop_upper(:,3);
    young_original = soil_prop_upper(:,2);
    
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
    
    k1=-10:20/1024:10;
    k2=-10:20/1024:10;
    
    nWave=length(k1);
    f = 0;
    nFreq=length(f);
    
    DAF_ratio = [10];
    
    zeroNod=find(feNod(:,3)==0); % find the index of the bottom nodes
    nEle = (length(zeroNod)-1)/2;
    nitmnode = zeros(nEle,1);
    
    for k=1:nEle
        nitmnode(k,1) = feNod(zeroNod(2*k),2); % find the coordinates of the middle nodes at the bottom
        midDof((1+3*(k-1)):3*k) = [6*k-2, 6*k-1, 6*k]; % find the middle nodes DOFs
    end
    
    
    zeroNod=find(feNod(:,3)==0);
    u=zeros(nFeDof,nWave);
    duy=zeros(nFeDof,nWave);
    u2=zeros(nFeDof,1);
    traction=zeros(dofs*nEle,nWave);
    
    % set up y range
    y0 = 0;
    
    for b = 1:length(DAF_ratio)
    c = DAF_ratio(b);
    [strcat('The current number of DAF value =', num2str(b))]
    
    for q = 1:itmax       % loop through iterations 
        
        % Define initial strain, stiffness & damping parameters
        % *******************************************************

        if (q == 1)  
            soil_prop_upper = soil_prop_upper_original ;
            soil_prop_inf = soil_prop_inf_original;
            Enew(:,1) = soil_prop_upper(:,2); % New calculated Young's modulus for each iteration
            Dnew(:,1) = soil_prop_upper(:,4);  % Damping for each iteration
            strain_new(:,1) = ones(length(element_number(:,1)),1); % Calculated ratios used for each iteration (i.e., G/Gmax)

        else           
            soil_prop_upper(:,2) = Enew(:,q);   % new thin layer properties of Young's modulus for equivalent stiffness and strain calculation
            soil_prop_inf(:,1) = Enew(1,q);   % new thin layer property of homogeneous halfspace (i.e., the property of the bottom thin layer)
            soil_prop_upper(:,4) = Dnew(:,q); % updated thin layer properties of damping
            soil_prop_inf(:,3) = Dnew(1,q);      
        end
        
        % Calculate green's function in each direction
        [uxxk,uyxk,uzxk,strain1x,strain2x,strain3x,strain4x,strain5x,strain6x,stress1x,stress2x,stress3x,stress4x,stress5x,stress6x]=greenfunction_PAC_full(k1,k2,cinf,f,c,nodes1,elements1,soil_prop_upper,soil_prop_inf,1);
        [uxyk,uyyk,uzyk,strain1y,strain2y,strain3y,strain4y,strain5y,strain6y,stress1y,stress2y,stress3y,stress4y,stress5y,stress6y]=greenfunction_PAC_full(k1,k2,cinf,f,c,nodes1,elements1,soil_prop_upper,soil_prop_inf,2);
        [uxzk,uyzk,uzzk,strain1z,strain2z,strain3z,strain4z,strain5z,strain6z,stress1z,stress2z,stress3z,stress4z,stress5z,stress6z]=greenfunction_PAC_full(k1,k2,cinf,f,c,nodes1,elements1,soil_prop_upper,soil_prop_inf,3);
        
        % calculate the traction at the bottom of the track and then return the displacement and the derivatives of the displacements of all the nodes
        [traction,uf,duf] = traction(nWave,nFreq,nFeDof,k1,k2,f,c,uxxk,uyxk,uzxk,uxyk,uyyk,uzyk,uxzk,uyzk,uzzk,nitmnode,zeroNod,feNod,midDof,dofs);
        % calculate the strains and stresses of each element in the track
        [epsilon,sigma] = track_element(nFeDof,nFeElt,dofs,nFreq,Ezz_s,nu_s,Ezz_a,nu_a,Ezz_b,nu_b,Ezz_usl,nu_usl,Ezz_sl,nu_sl,Ezz_rp,nu_rp,E_r,nu_r,uf,duf);
        % calculate the octahedral strains in the finite elements of track
        for iFreq=1:nFreq
            for iElt = 1:nFeElt
                oct_eps(iElt,iFreq)=octaheral(epsilon(:,iElt,iFreq),alpha);
            end
        end
        
        % calculate the strains in the free field in the subgrade
        [Strain_1, xx]=field_PAC_test(strain1z,strain1y,strain1x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Strain_2]=field_PAC_test(strain2z,strain2y,strain2x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Strain_3]=field_PAC_test(strain3z,strain3y,strain3x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Strain_4]=field_PAC_test(strain4z,strain4y,strain4x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Strain_5]=field_PAC_test(strain5z,strain5y,strain5x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Strain_6]=field_PAC_test(strain6z,strain6y,strain6x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        % calcuate the stresses in the free field in the subgrade
        [Stress_1]=field_PAC_test(stress1z,stress1y,stress1x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Stress_2]=field_PAC_test(stress2z,stress2y,stress2x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Stress_3]=field_PAC_test(stress3z,stress3y,stress3x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Stress_4]=field_PAC_test(stress4z,stress4y,stress4x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Stress_5]=field_PAC_test(stress5z,stress5y,stress5x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
        [Stress_6]=field_PAC_test(stress6z,stress6y,stress6x,k2,k1,nitmnode,zeroNod,feNod,y0,traction(3:3:end,:),traction(2:3:end-1,:),traction(1:3:end-2,:));
    
        confining_stress = track_weight/1000 + sigma1;  % confining stress = track weight + soil weight
        x_confi = strain_depth ;  % For plotting purpose
        oct = zeros(length(strain.Strain_1(:,1,1)),length(strain.Strain_1(1,1,:)));
        for w = 1:length(strain.Strain_1(1,1,:))
            oct(:,w) = alpha*sqrt((strain.Strain_1(:,1,w) - strain.Strain_2(:,1,w)).^2 + (strain.Strain_1(:,1,w) - strain.Strain_3(:,1,w)).^2 + (strain.Strain_2(:,1,w) - strain.Strain_3(:,1,w)).^2 + 6*(strain.Strain_4(:,1,w).^2 +strain.Strain_5(:,1,w).^2 +strain.Strain_6(:,1,w).^2 ) )./3;
            octmax(w) = max(oct(:,w));
            strain_val(q,w) = octmax(w);  % record the maximum (absolute value) of oactahedral shear strain in each thin layer
        end
        
        strain_val_upper1(q,:) = oct(:,5);
        
        [Gnew,Gratio,damping_modi] = modulus_reduction(young_original,strain_val(q,:),Ip,confining_stress);  % calculate new G (i.e., Young's modulus) and ratio (G/Gmax)
        
        Enew(:,q+1) = Gnew;  % update newly calculated Young's modulus
        strain_new(:,q+1) = Gratio;   % update newly calculated ratio (G/Gmax)
        Dnew(:,q+1) = damping_modi;  % updated damping info
        damp_ratio(q) = damping_modi(end);
        x_strain(q) = strain_val(q,end);  % record the strain/Gratio/damping values for plotting
        y_ratio(q) = Gratio(end);
        
        if ((abs((Enew(:,q+1) - Enew(:,q))./Enew(:,q)))<tol) == 1  % Program will end unless each of the thin layer satisfies the condition
            X = sprintf('The total number of iteration is %d', q);
            disp(X)     % print off the iteration results
            % Store the values for DAF curve and time history
            strain_DAF_upper1(b) = max(strain_val_upper1(end,:));
            
            clear uxxk uyxk uzxk strain1x strain2x strain3x strain4x strain5x strain6x stress1x stress2x stress3x stress4x stress5x stress6x
            clear uxyk uyyk uzyk strain1y strain2y strain3y strain4y strain5y strain6y stress1y stress2y stress3y stress4y stress5y stress6y
            clear uxzk uyzk uzzk strain1z strain2z strain3z strain4z strain5z strain6z stress1z stress2z stress3z stress4z stress5z stress6z
            clear Strain_1 Strain_2 Strain_3 Strain_4 Strain_5 Strain_6 Stress_1 Stress_2 Stress_3 Stress_4 Stress_5 Stress_6
            clear strain_new strain_val Enew Dnew damp_ratio x_strain y_ratio
            clear strain_val_upper1 

            break
        end
    end
    end
            
        
    xx1 = csvread('data/xx.csv');
    yy1 = csvread('data/yy.csv');
    zz1 = csvread('data/zz.csv');
    zx1 = csvread('data/zx.csv');

    % Note that the stresses at (0,0,-2) is compared, so we need to find the thin layer that corresponds to 2m 
    figure
    plot(xx,Stress_1(:,1,27)/1000,'r',xx,Stress_2(:,1,27)/1000,'b',xx,Stress_3(:,1,27)/1000,'g',xx,Stress_6(:,1,27)/1000,'m')
    hold on
    plot(xx1(:,1),-xx1(:,2),'*r')
    hold on
    plot(yy1(:,1),-yy1(:,2),'^b')
    hold on
    plot(zz1(:,1),-zz1(:,2),'+g')
    hold on
    plot(-zx1(:,1),-zx1(:,2),'om')
    axis([-25 25 -7 2])
    legend('TLM \sigma_{xx}','TLM \sigma_{yy}','TLM \sigma_{zz}','TLM \sigma_{zx}','Chen et al 2005 \sigma_{xx}','Chen et al 2005 \sigma_{yy}','Chen et al 2005 \sigma_{zz}','Chen et al 2005 \sigma_{zx}')
    xlabel('Distance along the track (m)')
    ylabel('Stress (kPa)')
    
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
    plot(1:20,traction(1:3:end-2,1),'b',1:20,traction(2:3:end-1,1),'r')
    
end  
    
    
    
    
    
    
    
    
    
    
    