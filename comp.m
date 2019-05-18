clear all
clc
close all

% addpath '.\asphalt_track_toolbox_v12'

dofs = 3; % 1: 1 dof; 2: 2 dofs; 3: 3 dofs

if(1) %computation?
    %% PROPERTIES
    
    %TRACK DATA
    %--------------------
    width_t=1.4/2;                                %<-- track: gauge [m]
    
    bl_s=5;                                         %<-- substrate: width (bottom) [m]
    h_s=0.4;                                       %<-- substrate: height [m]
    % slope_s=30;
    % bu_s=bl_s-2*h_s/tand(slope_s);                  %<-- substrate: width (upper) [m]
    bu_s=5;
    z0_s=0;                                         %<-- substrate: z coordinate at bottom [m]
    n1_s=50;%12;                                        %<-- substrate: number horizontal elements (mesh) [pair]
    n2_s=4;                                         %<-- substrate: number vertical elements (mesh)
    Ezz_s=200E6;                                    %<-- substrate: Young's modulus [N/m^2]
    nu_s=0.35;                                       %<-- substrate: Poisson's ratio [-]
    rho_s=2000;                                     %<-- substrate: density [kg/m^3]
    D_s=0.04;                                          %<-- substrate: damping ratio
    Ezz_s=Ezz_s*(1+sqrt(-1)*D_s);                   %<-- substrate: complex stiffness [N/m^2]
    
    
    bl_a=5;                                         %<-- substrate: width (bottom) [m]
    bu_a=5;
    n1_a=round(n1_s*bl_a/bu_s);                    %<-- asphalt: number horizontal elements (mesh) [pair]
    if(mod(n1_a,2));n1_a=n1_a+1;end
    bl_a=bu_s/n1_s*n1_a;
    bu_a=bl_a;
    Ezz_a=200E6;                                      %<-- asphalt: Young's modulus [N/m^2]
    nu_a=0.35;                                      %<-- asphalt: Poisson's ratio [-]
    rho_a=2000;                                     %<-- asphalt: density [kg/m^3]
    D_a=0.04*2;                                       %<-- asphalt: damping ratio
    Ezz_a=Ezz_a*(1+sqrt(-1)*D_a);                   %<-- asphalt: complex stiffness [N/m^2]
    h_a=0.1;                                       %<-- asphalt: height [m]
    n2_a=1;                                         %<-- substrate: number vertical elements (mesh)
    
    h_b=0.10;                                       %<-- ballast: height [m]
    n2_b=1;                                         %<-- ballast: number vertical elements (mesh)
    % slope_b=30;
    % bl_b=bu_s;                                      %<-- ballast: width (bottom) [m]
    % bu_b=bl_b-2*h_b/tand(slope_b);                  %<-- ballast: width (upper) [m]
    bl_b=5;                                      %<-- ballast: width (bottom) [m]
    
    n1_b=round(n1_a*bl_b/bu_a);
    if(mod(n1_b,2));n1_b=n1_b+1;end
    bl_b=bu_a/n1_a*n1_b;
    bu_b=bl_b;
    
    Ezz_b=200E6;                                    %<-- ballast: Young's modulus [N/m^2]
    nu_b=0.35;                                       %<-- ballast: Poisson's ratio [-]
    rho_b=2000;                                     %<-- ballast: density [kg/m^3]
    D_b=0.04;                                          %<-- ballast: damping ratio
    Ezz_b=Ezz_b*(1+sqrt(-1)*D_b);                   %<-- ballast: complex stiffness [N/m^2]
    
    d_sl=0.65;                                         %<-- sleeper: distance [m]
    Ezz_sl=30E9;                                      %<-- sleeper: Young's modulus [N/m^2]
    b_sl=0.25;
    l_sl=3.4;
    D_sl=0.06*2;
    Ezz_sl=Ezz_sl*(1+sqrt(-1)*D_sl);
    nu_sl=0.2;                                        %<-- sleeper: Poisson's ratio [-]
    h_sl=0.3;                                       %<-- sleeper: height [m]
    rho_sl=2400;                                       %<-- sleeper: height [m]
    n2_sl=3;                                         %<-- sleeper: number vertical elements (mesh)
    
    n1_sl=round(n1_b*l_sl/bu_b);
    if(mod(n1_sl,2));n1_sl=n1_sl+1;end
    l_sl=bu_b/n1_b*n1_sl;
    
    Ezz_usl=30E9;                                     %<-- under sleeper pad: Young's modulus [N/m^2]
    nu_usl=0.2;                                        %<--  under sleeper pad: Poisson's ratio [-]
    rho_usl=2400;                                       %<--  under sleeper pad: density [kg/m^3]
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
    
    kzz_rp=133E6;                                    %[N/m]
    b_rp=0.25;
    l_rp=0.1;
    h_rp=10e-3;                                       %<-- rail pad: height [m]
    
    Ezz_rp=8.889e6; %kzz_rp*h_rp/l_rp/b_rp;                                   %<-- rail pad: Young's modulus [N/m^2]
    nu_rp=0.49;                                        %<--  rail pad: Poisson's ratio [-]
    rho_rp=1000;                                       %<--  rail pad: density [kg/m^3]
    
    D_rp=0.021*2;                                    %<-- rail pad: damping [Ns^2/m]
    Ezz_rp=Ezz_rp*(1+sqrt(-1)*D_rp);                   %<-- rail pad: complex stiffness [N/m^2]
    
    Ezz_rp=Ezz_rp*1;                                   %<-- rail pad: stiffness [N/m/m]
    kzz_rp=kzz_rp*1;                                   %<-- rail pad: stiffness [N/m/m]
    
    h_r=0.155;
    l_r=l_rp;                
    A_r=0.0155;                                    %<-- rail: Section [m^2]
    Ixx_r=3.103e-5;                                %<-- rail: bending stiffness [m^4]-Ix
    Izz_r=1.292e-5;                                 % Iy
    E_r=2.1e11;                                     %<-- rail: Young's modulus [N/m^2]
    nu_r=0.3;                                       %<--  rail: Poisson's ratio [-]
    rho_r=7850;                                     %<-- rail: Density [kg/m^3]
    D_r=0.01;                                       %<-- rail: damping ratio
    E_r=E_r*(1+sqrt(-1)*D_r);                       %<-- rail: complex stiffness [N/m^2]
    rail_type = 1;                                  % rail type: 1-solid element; 2-beam element
    
    % SOIL DATA
    %------------------------
    % Soil density
    var1_d = {'p'};
    var1_v = [2000]';
    % Secondary wave speed
    var2_d = {'EE'};
    var2_v = [1e30];
    % Primary wave speed
    var3_d = {'v'};
    var3_v = [0.35]';
    % Soil layer heights
    layer_heights =[ 0.1]';
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
     
    
    %
    % fpp_emp = real(E_r*Ixx_r * (10.2*(sleeper_spacing^-1.61)*(m_sl^-0.33))^(1/0.33))    % pinned-pinned frequency
    % fs_emp = real(Ezz_b * (0.253*(sleeper_spacing^-0.15)*((E_r*Ixx_r)^-0.2))^2)          % Sleeper & rail resonance frequency
    % fr_emp = real(kzz_rp *(0.0275*(60^-0.3)*(m_sl^-0.11))^(1/0.59))    % Pad frequency
    %
    %
    % EI =(fpp_emp / (10.2*(L^-1.61)*(m^-0.33)))^(1/0.33);    % Rail EI
    % kb =(fs_emp / (0.253*(L^-0.15)*(EI^-0.2)))^2;           % Ballast stiffness
    % kp =(fr_emp /(0.0275*(m^-0.3)*(ms^-0.11)))^(1/0.59);    % Pad stiffness
    
    
    [feNod,feElt,feTyp,feMat,feSec,T1,t1,dofmaster,P]=asphalt_track_8nodes_com(width_t,bl_s,h_s,bu_s,z0_s,n1_s,n2_s,h_b,bu_b,bl_b,n2_b,n1_b,n2_a,n1_a,h_a,rho_s,Ezz_s,nu_s,...
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

    P1=T1*P;
    
    % SAMPLING
    % frequency sampling
    frequency=1:5:300;                           %<-- Receivers: Frequency sampling
    nFreq=length(frequency);
    omega=2*pi*frequency;
    % wavenumber sampling
    kdy=linspace(0,3.6,360);
    kdy=[kdy logspace(0.57415,2,32)];
    nWave=length(kdy);
    Cref=500; % the minumum secondary wave (Cs) speed in the soil
    py=kdy/Cref;
    c = 0;
    zeroNod=find(feNod(:,3)==0);
    
    %SOIL-STRUCTURE INTERACTION-PROBLEM
%     u=zeros(nFeDof+3,nWave,nFreq);
%     duy=zeros(nFeDof+3,nWave,nFreq);
    u=zeros(nFeDof,nWave);
    duy=zeros(nFeDof,nWave);
    unew = zeros(nFeDof,nWave,nFreq);
    duynew=zeros(nFeDof,nWave,nFreq);
    u2=zeros(nFeDof,1);
    
    
    tic;
    for iFreq=1:nFreq
%         iFreq
        u=zeros(nFeDof,nWave);
        duy=zeros(nFeDof,nWave);
        
        [uzk] = greenfunction(py.*omega(iFreq),cinf,c,nodes1,elements1,soil_prop_upper,soil_prop_inf); % calculate the greens function
        Ks = stiffvk(py.*omega(iFreq),width_t*2,uzk); % calculate the equivalent stiffness by eqn(5) in Costa's paper
           
        for iWave=1:nWave
            if omega(iFreq)==0
                ky=py(iWave);
            else
                ky=py(iWave)*omega(iFreq);
            end
            %FEM Matrices
            Kt=K0 + sqrt(-1)*omega(iFreq)*C0 - sqrt(-1)*K1*ky - K2*ky^2 + sqrt(-1)*K3*ky^3 + K4*ky^4 - omega(iFreq)^2*M;
            Kt=T1*Kt+t1;
            
            %Soil stiffness
            Kts=zeros(size(Kt,1),size(Kt,2));
            Kts(1:size(Kt,1),1:size(Kt,2))=Kt;
            
            switch dofs
                case 1
                    v = Ks(nWave)*ones(1,numel(zeroNod));
                    KtEq = zeros(size(Kt,1),size(Kt,2));
                    KtEq(1:numel*dofs,1:numel*dofs) = diag(v);
                    Kts = Kts+KtEq;
                case 2
                    A = [0 0; 0 Ks(nWave)];                                         % Original Matrix (Created)
                    N = numel(zeroNod);                                                  % Number Of Times To Repeat
                    Ar = repmat(A, 1, N);                                   % Repeat Matrix
                    Ac = mat2cell(Ar, size(A,1), repmat(size(A,2),1,N));    % Create Cell Array Of Orignal Repeated Matrix
                    Out = blkdiag(Ac{:});                                    % Desired Result
                    KtEq = zeros(size(Kt,1),size(Kt,2));
                    KtEq(1:numel(zeroNod)*dofs,1:numel(zeroNod)*dofs) = Out;
                    Kts = Kts+KtEq;
                case 3
                    A = [0 0 0; 0 0 0;0 0 Ks(nWave)];                                         % Original Matrix (Created)
                    N = numel(zeroNod);                                                  % Number Of Times To Repeat
                    Ar = repmat(A, 1, N);                                   % Repeat Matrix
                    Ac = mat2cell(Ar, size(A,1), repmat(size(A,2),1,N));    % Create Cell Array Of Orignal Repeated Matrix
                    Out = blkdiag(Ac{:});                                    % Desired Result
                    KtEq = zeros(size(Kt,1),size(Kt,2));
                    KtEq(1:numel(zeroNod)*dofs,1:numel(zeroNod)*dofs) = Out;
                    Kts = Kts+KtEq;
            end
            
            Kts = sparse(Kts);
            
            u2=Kts\P1;
%             u(1:1:end-3,iWave,iFreq)=u2;% depends on dofs
            u(1:1:end,iWave) = u2;
%             duy(1:1:end-3,iWave,iFreq)=-sqrt(-1)*ky*u2;
            duy(1:1:end,iWave)=-sqrt(-1)*ky*u2;
        end
        
        unew(:,:,iFreq) = u;% dispalcement in all dofs (include rail, railpad, etc etc) all freqs all wavenumbers
        duynew(:,:,iFreq) = duy; % duy = du/dy - derivative
    end
    time_tsi=toc;
    
    uf=zeros(nFeDof,nFreq);
    duyf=zeros(nFeDof,nFreq);
    duyf2=zeros(nFeDof,nFreq);
    
    %FREQUENCY RESPONSE
    y=0.001;
    for iDof=1:nFeDof
        for iFreq=1:nFreq
            if omega(iFreq)==0
                ky=py;
            else
                ky=py*omega(iFreq);
            end
            uf(iDof,iFreq)=1/pi*intfilon(ky,-sqrt(-1)*y,unew(iDof,:,iFreq)); % transform to frequency domain only
            duf(iDof,iFreq)=1/pi*intfilon(ky,-sqrt(-1)*y,duynew(iDof,:,iFreq));
        end
    end
    
    figure
    plot(frequency,abs(uf(end-4*dofs,:))) % plot displacement on rail only
    hold on
    plot(frequency,abs(uf(end-dofs,:))) 
    
    
    

vidObj = VideoWriter('disp_sab_v5_DC.avi');
vidObj.Quality=100;
open(vidObj);
figure1 = figure('Renderer','zbuffer','Color',[1 1 1]);
% ax1 = axes('Parent',figure1,'Visible','off','Fontsize',20);
warning off
for iFreq=1:nFreq
    X=feNod(:,2);
    Y=feNod(:,3);
    tri=delaunay(X,Y);
    trisurf(tri,X,Y,abs(squeeze(uf(3:3:end,iFreq))));
    shading('interp')
    view([0 -90]);
    hold on
    plot(X,Y,'ok')
    set (gca,'Ydir','reverse')
    hold off
    xlabel ('X [m]')
    ylabel ('Y [m]')
    text(0.9,0.8,'Ballast','Units','normalized')
    text(0.9,0.4,'Asphalt','Units','normalized')
    text(0.9,0.15,'Subgrade','Units','normalized')
    %     legend(strcat(int2str(frequency(iFreq)),'Hz'));
    %     title(strcat('Displacement [m/kN] Frequency ',int2str(frequency(iFreq)),'Hz'));
%     ax1.Visible = 'off';
%     ax1.XTick = [];
%     ax1.YTick = [];
    % colormap
%     colormap(ax1,'jet')
    % colorbar
    c = colorbar;
%     caxis([min(abs(squeeze(uf(3:3:end,iFreq)))) max(abs(squeeze(uf(3:3:end,iFreq))))])
    drawnow
    F = getframe(gcf);
    writeVideo(vidObj,F);
end
warning on
close(vidObj);

recep_dc_freq = recep_02(:,1);
recep_dc_mag = recep_02(:,2);

figure
plot(recep_dc_freq,recep_dc_mag)
xlim([0 300])
hold on
plot(recep_kd_freq,recep_kd_mag)
xlabel('Frequency Hz')
ylabel('Receptance')
legend('DC','KD')

recep_kd_freq = frequency;
recep_kd_mag = abs([3.83185970475275e-10 - 1.87431325713781e-11i,8.96311322987992e-10 - 4.01135299394943e-11i,9.23715634317490e-10 - 4.05816666496068e-11i,9.34316497532822e-10 - 4.09567020318463e-11i,9.42964363044217e-10 - 4.14175597084742e-11i,9.51562194596082e-10 - 4.19979938924632e-11i,9.60680451024245e-10 - 4.27235577885191e-11i,9.70630303941682e-10 - 4.36236702756570e-11i,9.81703514193782e-10 - 4.47368215393824e-11i,9.94255088702130e-10 - 4.61272909906809e-11i,1.00873296362897e-09 - 4.79016956598490e-11i,1.02573916185934e-09 - 5.02451917656731e-11i,1.04616618799485e-09 - 5.35100241267380e-11i,1.07147910203825e-09 - 5.84717505885861e-11i,1.10430078358723e-09 - 6.73162378705561e-11i,1.14791322273050e-09 - 8.58639143680642e-11i,1.20346528955792e-09 - 1.36477655371610e-10i,1.20819383464435e-09 - 2.42639876240084e-10i,1.13024631759850e-09 - 2.83391417417030e-10i,1.09656058078313e-09 - 2.51806631225662e-10i,1.11340559036488e-09 - 2.17667664435230e-10i,1.16211974149166e-09 - 1.96143872517678e-10i,1.23976265161452e-09 - 1.89544671600405e-10i,1.35669021124160e-09 - 2.03634485500978e-10i,1.54417472820520e-09 - 2.60621450808171e-10i,1.88761797054319e-09 - 4.70109009147809e-10i,2.15719333293741e-09 - 1.65378639050106e-09i,6.87923477632958e-10 - 1.42180750882993e-09i,3.86583862537751e-10 - 1.68966487949588e-09i,-3.04059870741542e-11 - 1.04902488452740e-09i,1.40705063874044e-10 - 6.20730935379516e-10i,2.92733196476165e-10 - 4.66237108842076e-10i,3.92036541604444e-10 - 4.06799202642875e-10i,4.52510153415352e-10 - 3.74715474298902e-10i,4.93318703029715e-10 - 3.57220495026919e-10i,5.17417348397317e-10 - 3.45723300114969e-10i,5.31539485052569e-10 - 3.33133498431478e-10i,5.41522411811927e-10 - 3.18487978991202e-10i,5.50358947521392e-10 - 3.02308947160003e-10i,5.59437873752417e-10 - 2.85659604050516e-10i,5.69018758803237e-10 - 2.68695673768544e-10i,5.79737755960189e-10 - 2.50991408856661e-10i,5.94147411124017e-10 - 2.31116657099711e-10i,6.19059602313980e-10 - 2.14477276325120e-10i,6.43405944346009e-10 - 2.16100577134204e-10i,6.53667181705285e-10 - 2.16113034294099e-10i,6.63773412619130e-10 - 2.11174802179805e-10i,6.76058973380188e-10 - 2.04388394020627e-10i,6.91165995441041e-10 - 1.98725556609275e-10i,7.07175226593270e-10 - 1.95145454417713e-10i,7.23904885289887e-10 - 1.92710490337072e-10i,7.42299584785434e-10 - 1.93334367621824e-10i,7.58535154216821e-10 - 1.96371388300686e-10i,7.76857326056936e-10 - 2.00667530072906e-10i,7.93517551150967e-10 - 2.14060672804681e-10i,7.93022445468470e-10 - 2.27864185262285e-10i,7.92702545571037e-10 - 2.23228642457675e-10i,8.03360142995373e-10 - 2.16069227532716e-10i,8.19674908573218e-10 - 2.11681768062274e-10i,8.40125082379669e-10 - 2.11572358818249e-10i]);

recep_kd_freq02 = frequency;
recep_kd_mag02 = abs([9.39149340169271e-11 - 4.17261305048851e-12i,1.42739198828428e-10 - 7.87712334988158e-12i,1.43544868636398e-10 - 7.97146387535709e-12i,1.44843426400098e-10 - 8.12103922220708e-12i,1.46670437394681e-10 - 8.33163053841429e-12i,1.49069927674299e-10 - 8.61180681693300e-12i,1.52090137126105e-10 - 8.97353653282059e-12i,1.55792096326914e-10 - 9.43275795642075e-12i,1.60277966732350e-10 - 1.00115724819842e-11i,1.65711016702918e-10 - 1.07470522082841e-11i,1.72325620339245e-10 - 1.17012779391315e-11i,1.80460369201949e-10 - 1.29829233180991e-11i,1.90639277149772e-10 - 1.48010527951662e-11i,2.03752688859188e-10 - 1.76275965462995e-11i,2.21461976671713e-10 - 2.29992653632517e-11i,2.44540222326625e-10 - 3.37821752774202e-11i,2.75871614560270e-10 - 6.30697761397088e-11i,2.77824345098720e-10 - 1.26119003805191e-10i,2.30096719590216e-10 - 1.52243659140123e-10i,2.06515513124212e-10 - 1.35723009182659e-10i,2.11531623095085e-10 - 1.16888744386133e-10i,2.34430473257044e-10 - 1.04732725949659e-10i,2.73233916666572e-10 - 1.00605438423884e-10i,3.33279408925511e-10 - 1.07571905391137e-10i,4.31054632963727e-10 - 1.37333880861157e-10i,6.11623366998736e-10 - 2.47616296804012e-10i,7.54079749731228e-10 - 8.69003803875657e-10i,3.27864839597228e-12 - 7.50601549995551e-10i,3.80008492726785e-11 - 9.57613617857601e-10i,-4.41291055444194e-10 - 5.70459275064174e-10i,-3.30027622690852e-10 - 3.33137353511038e-10i,-2.35769637077224e-10 - 2.67685310489131e-10i,-2.02046511966375e-10 - 2.64976585334652e-10i,-1.99096551492395e-10 - 2.46109716929604e-10i,-1.94828186148094e-10 - 2.26113063457574e-10i,-1.92055608697726e-10 - 2.08499749959215e-10i,-1.90254554468170e-10 - 1.91955654454286e-10i,-1.88457476501977e-10 - 1.76256582160132e-10i,-1.86181917036574e-10 - 1.61469009051875e-10i,-1.83192297098962e-10 - 1.47715704372280e-10i,-1.79298997136771e-10 - 1.35059533806100e-10i,-1.74183366201598e-10 - 1.23692131408054e-10i,-1.67018654145062e-10 - 1.14060420544186e-10i,-1.56366840554716e-10 - 1.10258643460377e-10i,-1.55639407345618e-10 - 1.18015595605820e-10i,-1.64428816395563e-10 - 1.12900705631490e-10i,-1.66476948720980e-10 - 1.03626748042896e-10i,-1.64715980042488e-10 - 9.62592300147107e-11i,-1.63163498450102e-10 - 9.16271491171876e-11i,-1.61416414191984e-10 - 8.76767986185455e-11i,-1.58848430997095e-10 - 8.54514236976140e-11i,-1.57047547691353e-10 - 8.95967113464561e-11i,-1.67254005197999e-10 - 9.30295365692428e-11i,-1.75536937237781e-10 - 8.62746529977254e-11i,-1.78215432841586e-10 - 8.03746922445542e-11i,-1.81538000409145e-10 - 7.86856337729588e-11i,-1.88849549880449e-10 - 7.65544436553231e-11i,-1.97686978949141e-10 - 7.17631244248739e-11i,-2.06820625117557e-10 - 6.40023647086012e-11i,-2.14409903638617e-10 - 5.23903781077062e-11i]);

recep_dc_freq2 = recep(:,1);
recep_dc_mag2 = recep(:,2);

figure
plot(recep_dc_freq2,recep_dc_mag2)
xlim([0 300])
hold on
plot(recep_kd_freq02,recep_kd_mag02)
xlabel('Frequency Hz')
ylabel('Receptance')
legend('DC','KD')

    figure
    X=feNod(:,2);
    Y=feNod(:,3);
    tri=delaunay(X,Y);
    trisurf(tri,X,Y,abs(uf(3:3:end,60)));
    shading('interp')
    view([0 -90]);
    hold on
    plot(X,Y,'ok')
    set (gca,'Ydir','reverse')
    hold off
    xlabel ('X [m]')
    ylabel ('Y [m]')
    c = colorbar;
    drawnow
    
    %STRAIN STRESS SOLID ELEMENT RESULTS
    
    % 6 strain/stress elements: xx,yy,zz,xy,zy,zx
    
    stress=1;
    sigman=zeros(6,nFeDof/dofs,nFreq); % depends on dofs
    sigman2=zeros(6,nFeElt,nFeDof/dofs,nFreq);
    epsilon=zeros(6,8,nFeElt,nFreq);% first value 6 means the number of strain elements                                      %[epsilon_zz; epsilon_xz; epsilonyz]
    sigma=zeros(6,8,nFeElt,nFreq); 
    
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
                    [epsilon(:,:,iElt,iFreq),sigma(:,:,iElt,iFreq)]=volu8s(nodcoord,uf(nodes,iFreq),duf(nodes,iFreq),E,nu); % epsilon xy,xz,yz shall be divided by 2 for the real values
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
                    sigman2(:,iElt,nodes,iFreq)=squeeze(sigma(:,:,iElt,iFreq));
                elseif(min(feElt(iElt,5:12))==0)
                    nodes=feElt(iElt,5:11);
                    sigman2(:,iElt,nodes,iFreq)=squeeze(sigma(:,1:7,iElt,iFreq));
                end
            end
        end
        
        % in case the node that is shared by multiple elements have
        % different stresses, the code below provides a solution to the
        % situation by taking the average of the stresses calculated by
        % those elements
        
        % in most cases, sigman2 = sigman
        for istress=1:6 
            for iFreq=1:nFreq
                for inod=1:1:nFeDof/dofs
                    nm=length(nonzeros(sigman2(istress,:,inod,iFreq)));
                    sigman(istress,inod,iFreq)=sum(nonzeros(sigman2(istress,:,inod,iFreq)))/nm;
                end
            end
        end
        clear sigman2
    end
    
    
    % elements
    iElts=75;
  
    
    istress1=1;
    istress2=2;
    istress3=3;
    inode=3;
    
    figure
    hold on
    plot(frequency,abs(squeeze(sigma(istress1,inode,iElts,:))),'-b')
    plot(frequency,abs(squeeze(sigma(istress2,inode,iElts,:))),'-r')
    plot(frequency,abs(squeeze(sigma(istress3,inode,iElts,:))),'-k')
    plot(frequency,abs(squeeze(sigma(4,inode,iElts,:))),'--')
    plot(frequency,abs(squeeze(sigma(5,inode,iElts,:))),'--')
    plot(frequency,abs(squeeze(sigma(6,inode,iElts,:))),'--')
    xlabel('Frequency (Hz)')
    ylabel('Stress [N/m^2/N]')
    legend('xx','yy','zz','xy','xz','yz')
    
    % Principal stresses
    sigma_mat=zeros(3,3,nFreq);
    sigma1=zeros(nFreq,1);sigma2=sigma1;sigma3=sigma1;
    for iFreq=1:nFreq
        
        sigma_mat(:,:,iFreq) = [sigma(1,inode,iElts,iFreq),sigma(4,inode,iElts,iFreq),sigma(5,inode,iElts,iFreq);
            sigma(4,inode,iElts,iFreq),sigma(2,inode,iElts,iFreq),sigma(6,inode,iElts,iFreq);
            sigma(5,inode,iElts,iFreq),sigma(6,inode,iElts,iFreq),sigma(3,inode,iElts,iFreq)];
   
        d=eigs(sigma_mat(:,:,iFreq));
        sigma1(iFreq,1) = d(1);
        sigma2(iFreq,1) = d(2);
        sigma3(iFreq,1) = d(3);
    end
    
    figure
    plot(frequency,abs(sigma1),'r',frequency,abs(sigma2),'b',frequency,abs(sigma3),'g')
    legend('S11','S22','S33')
    
    istress=1;
    figure
    X=feNod(:,2);
    Y=feNod(:,3);
    tri=delaunay(X,Y);
    trisurf(tri,X,Y,abs(sigman(istress,:,1)));
    shading('interp')
    view([0 -90]);
    hold on
    plot(X,Y,'ok')
    set (gca,'Ydir','reverse')
    hold off
    xlabel ('X [m]')
    ylabel ('Y [m]')
    c = colorbar;
    drawnow
    
        
      
    
    
    
    
    
    
end 
    