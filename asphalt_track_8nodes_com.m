%--------------------------------------------------------------
% asphalt_track.m returns the node,element, force information in the track
% model
%--------------------------------------------------------------
% Output:
%
% feNod: Complete nodes information matrix (size = n_node * 3, n_node means the total number of nodes)
%        1st column: Node number
%        2nd column: x-coordinate of the corresponding node
%        3rd column: y-coordinate of the corresponding node
%
% feElt: Complete element information matrix (size = n_ele * 8, n_ele means the total number of elements)
%        1st column: Element number
%        2nd-4th column: Label of the structure (i.e., 1-subballast, 2 for asphalt, etc)
%        5th-8th column: The corresponding nodes number in the element, going anti-clockwise
%
% feTyp: Type of the structures (size = n_struct * 2)
%        1st column: Number of the structures, which should match with the number in the 2nd-4th column of feElt
%        2nd column: Type of the structure - e.g 'volu4' or 'beam1' etc
%
% feMat: Properties matrix of each structure (size = n_struct * 4)
%        1st column: Number of the structures, which should match with the number in the 2nd-4th column of feElt
%        2nd column: Young's modulus
%        3rd column: Poisson's ratio
%        4th column: Density
%
% feSec: Geometric data matrix of the section (could be any size)
%        It can be used to store the area/mass/etc information of the section
%
% T1, t1: Constraint matrices
%
% nodemaster: The node number for the masternode
% 
% P: The position of the force acted on the rail nodes
%
% Model description: It contains 7 parts in this model, from bottom to top,
%                    subballast, asphalt, ballast, sleeper pad, sleeper,
%                    railpad, rail. The information matrix of each part is
%                    added together to form the global property matrix.


function [feNod,feElt,feTyp,feMat,feSec,T1,t1,nodemaster,P]=asphalt_track_8nodes_com(width_t,bl_s,h_s,bu_s,z0_s,n1_s,n2_s,h_b,bu_b,bl_b,n2_b,n1_b,n2_a,n1_a,h_a,rho_s,Ezz_s,nu_s,...
    rho_b,Ezz_b,nu_b,rho_a,Ezz_a,nu_a,bl_a,bu_a,rho_usl,Ezz_usl,nu_usl,h_usl,n2_usl,n1_usl,l_usl,rho_sl,Ezz_sl,nu_sl,h_sl,n2_sl,n1_sl,l_sl,...
    Ezz_rp,rho_rp,nu_rp,h_rp,l_rp,E_r,l_r,h_r,rho_r,nu_r,A_r,Ixx_r,Izz_r,dofs,rail_type)


% FINITE ELEMENT MODEL

% subgrade bottom
[feNod_s,feElt_s2]=solidtrack_8nodes_new(bl_s,bu_s,n1_s,h_s,n2_s,width_t,z0_s,l_rp,0);
nNod_s=length(feNod_s);
nElt_s=length(feElt_s2);
feElt_s=ones(nElt_s,12);
feElt_s(:,1)=feElt_s2(1:nElt_s,1);
feElt_s(:,5:12)=feElt_s2(1:nElt_s,2:9);
% sevde = feElt_s(:,end)==0;
% feElt_s(sevde,2)=2;

plot(feNod_s(:,2),feNod_s(:,3),'or')

%sleeper

z0_sl=h_s;
[feNod_sl2,feElt_sl2]=solidtrack_8nodes_new(l_sl,l_sl,n1_sl,h_sl,n2_sl,width_t,z0_sl,l_rp,0);
nNod_sl=length(feNod_sl2);
nElt_sl=length(feElt_sl2);
feNod_sl=feNod_sl2;
feNod_sl(:,1)=feNod_sl(:,1)+nNod_s;
feElt_sl=ones(nElt_sl,12);
feElt_sl(:,1)=feElt_sl2(1:nElt_sl,1)+nElt_s;
feElt_sl(:,5:12)=feElt_sl2(1:nElt_sl,2:9)+nNod_s;
feElt_sl(:,2:4)=2;
hold on
plot(feNod_sl(:,2),feNod_sl(:,3),'^m')

%rail pad
[~,p]=min(abs(feNod_sl(end-2*n1_sl:end,2)+width_t));
nNod_rp=16;
nElt_rp=2;
feNod_rp(1:nNod_rp,1)=1:1:nNod_rp;
feNod_rp(1:nNod_rp,1)=feNod_rp(1:nNod_rp,1)+nNod_s+nNod_sl;
feNod_rp(1:6,3)=h_s+h_sl;
feNod_rp(7:10,3)=h_s+h_sl+h_rp/2;
feNod_rp(11:16,3)=h_s+h_sl+h_rp;

feNod_rp(1:3,2)=[feNod_sl(end-2*n1_sl+p-1,2)-l_rp/2 feNod_sl(end-2*n1_sl+p-1,2) feNod_sl(end-2*n1_sl+p-1,2)+l_rp/2];
feNod_rp(4:6,2)=[feNod_sl(end-p+1,2)-l_rp/2 feNod_sl(end-p+1,2) feNod_sl(end-p+1,2)+l_rp/2];
feNod_rp(7:8,2)=[feNod_sl(end-2*n1_sl+p-1,2)-l_rp/2 feNod_sl(end-2*n1_sl+p-1,2)+l_rp/2];
feNod_rp(9:10,2)=[feNod_sl(end-p+1,2)-l_rp/2 feNod_sl(end-p+1,2)+l_rp/2];
feNod_rp(11:13,2)=[feNod_sl(end-2*n1_sl+p-1,2)-l_rp/2 feNod_sl(end-2*n1_sl+p-1,2) feNod_sl(end-2*n1_sl+p-1,2)+l_rp/2];
feNod_rp(14:16,2)=[feNod_sl(end-p+1,2)-l_rp/2 feNod_sl(end-p+1,2) feNod_sl(end-p+1,2)+l_rp/2];

% 1-5-2-6-3-7-4-8
feElt_rp=[1+nElt_s+nElt_sl 3 3 3 feNod_rp(1,1) feNod_rp(3,1) feNod_rp(13,1) feNod_rp(11,1) feNod_rp(2,1) feNod_rp(8,1) feNod_rp(12,1) feNod_rp(7,1);
    2+nElt_s+nElt_sl 3 3 3 feNod_rp(4,1) feNod_rp(6,1) feNod_rp(16,1) feNod_rp(14,1) feNod_rp(5,1) feNod_rp(10,1) feNod_rp(15,1) feNod_rp(9,1)];


% feElt_rp=[1+nElt_s+nElt_sl 3 3 3 feNod_rp(1,1) feNod_rp(2,1) feNod_rp(3,1) feNod_rp(8,1) feNod_rp(13,1) feNod_rp(12,1) feNod_rp(11,1) feNod_rp(7,1);
%     2+nElt_s+nElt_sl 3 3 3 feNod_rp(4,1) feNod_rp(5,1) feNod_rp(6,1) feNod_rp(10,1) feNod_rp(16,1) feNod_rp(15,1) feNod_rp(14,1) feNod_rp(9,1)];

hold on
plot(feNod_rp(:,2),feNod_rp(:,3),'*k')

%rail
if rail_type==1
    nNod_r = 26;
    feNod_r(1:nNod_r,1)=1:1:nNod_r;
    feNod_r(1:nNod_r,1)=feNod_r(1:nNod_r,1)+nNod_s+nNod_sl+nNod_rp;
    feNod_r(1:6,3)=feNod_rp(11:16,3);
    feNod_r(11:16,3)=feNod_rp(11:16,3)+h_r/2;
    feNod_r(7:10,3)=feNod_rp(7:10,3)+h_r/4+h_rp/2;
    feNod_r(17:20,3)=feNod_r(7:10,3)+h_r/2;
    feNod_r(21:26,3)=feNod_r(1:6,3)+h_r;
    
    feNod_r(1:6,2)=feNod_rp(11:16,2);
    feNod_r(11:16,2)=feNod_r(1:6,2);
    feNod_r(21:26,2)=feNod_r(11:16,2);
    feNod_r(7:10,2)=feNod_rp(7:10,2);
    feNod_r(17:20,2)=feNod_r(7:10,2);
   
    %  1-2-3-4-5-6-7-8
%     feElt_r=[1+nElt_s+nElt_sl+nElt_rp 4 4 4 feNod_r(1,1) feNod_r(2,1) feNod_r(3,1) feNod_r(8,1) feNod_r(13,1) feNod_r(12,1) feNod_r(11,1) feNod_r(7,1);
%     2+nElt_s+nElt_sl+nElt_rp 4 4 4 feNod_r(4,1) feNod_r(5,1) feNod_r(6,1) feNod_r(10,1) feNod_r(16,1) feNod_r(15,1) feNod_r(14,1) feNod_r(9,1);
%     3+nElt_s+nElt_sl+nElt_rp 4 4 4 feNod_r(11,1) feNod_r(12,1) feNod_r(13,1) feNod_r(18,1) feNod_r(23,1) feNod_r(22,1) feNod_r(21,1) feNod_r(17,1);
%     4+nElt_s+nElt_sl+nElt_rp 4 4 4 feNod_r(14,1) feNod_r(15,1) feNod_r(16,1) feNod_r(20,1) feNod_r(26,1) feNod_r(25,1) feNod_r(24,1) feNod_r(19,1);];

    % 1-5-2-6-3-7-4-8
    feElt_r=[1+nElt_s+nElt_sl+nElt_rp 4 4 4 feNod_r(1,1) feNod_r(3,1) feNod_r(13,1) feNod_r(11,1) feNod_r(2,1) feNod_r(8,1) feNod_r(12,1) feNod_r(7,1);
    2+nElt_s+nElt_sl+nElt_rp 4 4 4 feNod_r(4,1) feNod_r(6,1) feNod_r(16,1) feNod_r(14,1) feNod_r(5,1) feNod_r(10,1) feNod_r(15,1) feNod_r(9,1);
    3+nElt_s+nElt_sl+nElt_rp 4 4 4 feNod_r(11,1) feNod_r(13,1) feNod_r(23,1) feNod_r(21,1) feNod_r(12,1) feNod_r(18,1) feNod_r(22,1) feNod_r(17,1);
    4+nElt_s+nElt_sl+nElt_rp 4 4 4 feNod_r(14,1) feNod_r(16,1) feNod_r(26,1) feNod_r(24,1) feNod_r(15,1) feNod_r(20,1) feNod_r(25,1) feNod_r(19,1);];

    
    hold on
    plot(feNod_r(:,2),feNod_r(:,3),'^r')
elseif rail_type==2
    feNod_r(1,1)=1+nNod_s+nNod_a+nNod_b+nNod_usl+nNod_sl+nNod_rp;
    feNod_r(2,1)=2+nNod_s+nNod_a+nNod_b+nNod_usl+nNod_sl+nNod_rp;
    feNod_r(1,2)=-width_t; %feNod_rp(2,2); 
    feNod_r(2,2)=+width_t; %feNod_rp(5,2); 
    feNod_r(:,3)=h_s+h_a+h_b+h_usl+h_sl+h_rp;
    feElt_r=[1+nElt_s+nElt_b+nElt_a+nElt_usl+nElt_sl+nElt_rp 7 7 7 feNod_r(1,1) 0 0 0 0 0 0 0 ;
         2+nElt_s+nElt_b+nElt_a+nElt_usl+nElt_sl+nElt_rp 7 7 7 feNod_r(2,1) 0 0 0 0 0 0 0 ];
    hold on
    plot(feNod_r(:,2),feNod_r(:,3),'^r')
end

if rail_type==1 % solid element
    
    feTyp={1 'volu8';2 'volu8';3 'volu8';4 'volu8'};
    feSec=[];
    feMat=[1 Ezz_s nu_s rho_s;
    2 Ezz_sl nu_sl rho_sl;
    3 Ezz_rp nu_rp rho_rp;
    4 E_r nu_r rho_r];

elseif rail_type==2  % beam element
    
    feTyp={1 'volu8';2 'volu8';3 'volu8';4 'volu8'; 5 'volu8';6 'volu8';7 'volu8';8 'beam1'};
    feSec=[8 A_r Ixx_r Izz_r];
    feMat=[1 Ezz_s nu_s rho_s;
    2 Ezz_s nu_s rho_s;
    3 Ezz_a nu_a rho_a;
    4 Ezz_b nu_b rho_b;
    5 Ezz_usl nu_usl rho_usl;
    6 Ezz_sl nu_sl rho_sl;
    7 Ezz_rp nu_rp rho_rp;
    8 E_r 0 rho_r];

end
        


feElt=[feElt_s;feElt_sl;feElt_rp;feElt_r];
feNod=[feNod_s;feNod_sl;feNod_rp;feNod_r];

%LOAD ON RAIL
nFeDof=dofs*size(feNod,1);% depends on dof
P=zeros(nFeDof,1);
if rail_type==1
    P(end-dofs)=0.5;
    P(end-4*dofs)=0.5; % depends on the dof, (en3-2)*nof+nof

elseif rail_type==2
    P(end)=0.5;
    P(end-dofs)=0.5;
end
% plot(feNod(:,2),feNod(:,3),'or')

%coupled nodes are identified

%soil-track interface
cpind=[1:1:(2*n1_s+1)]';
indnm=find(feNod(cpind,2)==0);
nodemaster=feNod(indnm,1);
cpind(nodemaster)=[];
n_nod_coupled_ss=2*n1_s;

%subgrade-asphalt interface

n_nod_coupled_ssl=2*n1_sl+1;
coup_nod_ssl=zeros(n_nod_coupled_ssl,2);
coup_nod_ssl(:,1)=feNod_sl(1:2*n1_sl+1,1);
for inode=1:n_nod_coupled_ssl
    coup_nod_ssl(inode,2)=feNod_s(find((abs(feNod_s(:,2)-feNod(coup_nod_ssl(inode,1),2))<1e-10)&(abs(feNod_s(:,3)-feNod(coup_nod_ssl(inode,1),3))<1e-10)),1);
end

%rail pad-sleeper
n_nod_coupled_rps=6;
coup_nod_rps=zeros(n_nod_coupled_rps,2);
coup_nod_rps(:,1)=feNod_rp(1:6,1);
 for inode=1:n_nod_coupled_rps
        coup_nod_rps(inode,2)=feNod_sl(find((abs(feNod_sl(:,2)-feNod(coup_nod_rps(inode,1),2))<1e-10)&(abs(feNod_sl(:,3)-feNod(coup_nod_rps(inode,1),3))<1e-10)),1);
 end


% for inode=1:n_nod_coupled_rps/2
%     coup_nod_rps(inode,2)=feNod_sl(end-n1_sl+p-1,1);
% end
% 
% for inode=n_nod_coupled_rps/2+1:n_nod_coupled_rps
%     coup_nod_rps(inode,2)=feNod_sl(end-p+1,1);
% end

%rail-rail pad
if rail_type==1
    n_nod_coupled_rpr=6;
    
    coup_nod_rpr=zeros(n_nod_coupled_rps,2);
    coup_nod_rpr(:,1)=feNod_r(1:6,1);
    for inode=1:n_nod_coupled_rpr
        coup_nod_rpr(inode,2)=feNod_rp(find((abs(feNod_rp(:,2)-feNod(coup_nod_rpr(inode,1),2))<1e-10)&(abs(feNod_rp(:,3)-feNod(coup_nod_rpr(inode,1),3))<1e-10)),1);
    end
    
elseif rail_type==2
    coup_nod_rpr(1:3,1)=feNod_r(1,1);
    coup_nod_rpr(1,2)=feNod_rp(11,1);
    coup_nod_rpr(2,2)=feNod_rp(12,1);
    coup_nod_rpr(3,2)=feNod_rp(13,1);
    coup_nod_rpr(4:6,1)=feNod_r(2,1);
    coup_nod_rpr(4,2)=feNod_rp(14,1);
    coup_nod_rpr(5,2)=feNod_rp(15,1);
    coup_nod_rpr(6,2)=feNod_rp(16,1);
     
    n_nod_coupled_rpr=6;
end
%rail-rail
if rail_type==1
    n_nod_coupled_rr=4;
    coup_nod_rr=zeros(n_nod_coupled_rr,2);
    nodrail1=[11 13 14 16];
    nodrail2=[12 12 15 15];
    coup_nod_rr(:,1)=feNod_r(nodrail1,1);
    coup_nod_rr(:,2)=feNod_r(nodrail2,1);
end

%Add constraint equations
%soil-track interface
Constr_ss=zeros(n_nod_coupled_ss*dofs,5); % multiply by 3 means 3 dofs
Constr_ss(:,2)=-1;
Constr_ss(:,4)=1;
switch dofs
    case 1
        Constr_ss(:,3)=cpind;
        Constr_ss(:,5)=nodemaster;
    case 2
        cpind2=zeros(length(cpind)*2,1);
        for icpind=1:length(cpind)
            cpind2((icpind-1)*2+1)=(cpind(icpind)-1)*2+1;
            cpind2((icpind-1)*2+2)=(cpind(icpind)-1)*2+2;
        end
        Constr_ss(:,3)=cpind2;
        Constr_ss(:,5)=repmat([(nodemaster-1)*2+1; (nodemaster-1)*2+2],[n_nod_coupled_ss 1]);
    case 3
        cpind2=zeros(length(cpind)*3,1); % multiply by 3 means 3 dofs
        for icpind=1:length(cpind)
            cpind2((icpind-1)*3+1)=(cpind(icpind)-1)*3+1;
            cpind2((icpind-1)*3+2)=(cpind(icpind)-1)*3+2;% multiply by 3 means 3 dofs
            cpind2((icpind-1)*3+3)=(cpind(icpind)-1)*3+3;
        end
        Constr_ss(:,3)=cpind2;
        Constr_ss(:,5)=repmat([(nodemaster-1)*3+1; (nodemaster-1)*3+2; (nodemaster-1)*3+3],[n_nod_coupled_ss 1]);% multiply by 3 means 3 dofs
end
%subgrade-asphalt interface
Constr_ssl=zeros(n_nod_coupled_ssl*dofs,5); % multiply by 3 means 3 dofs 
Constr_ssl(:,2)=-1;
Constr_ssl(:,4)=1;
for i_constr=1:n_nod_coupled_ssl
    switch dofs
        case 1
            Constr_ssl((i_constr-1)*1+1,3)=coup_nod_ssl(i_constr,1);
            Constr_ssl((i_constr-1)*1+1,5)=coup_nod_ssl(i_constr,2);
        case 2
            Constr_ssl((i_constr-1)*2+1,3)=(coup_nod_ssl(i_constr,1)-1)*2+1;
            Constr_ssl((i_constr-1)*2+1,5)=(coup_nod_ssl(i_constr,2)-1)*2+1;
            Constr_ssl((i_constr-1)*2+2,3)=(coup_nod_ssl(i_constr,1)-1)*2+2;
            Constr_ssl((i_constr-1)*2+2,5)=(coup_nod_ssl(i_constr,2)-1)*2+2;
        case 3           
            Constr_ssl((i_constr-1)*3+1,3)=(coup_nod_ssl(i_constr,1)-1)*3+1;
            Constr_ssl((i_constr-1)*3+1,5)=(coup_nod_ssl(i_constr,2)-1)*3+1;% multiply by 3 means 3 dofs
            Constr_ssl((i_constr-1)*3+2,3)=(coup_nod_ssl(i_constr,1)-1)*3+2;
            Constr_ssl((i_constr-1)*3+2,5)=(coup_nod_ssl(i_constr,2)-1)*3+2; 
            Constr_ssl((i_constr-1)*3+3,3)=(coup_nod_ssl(i_constr,1)-1)*3+3;
            Constr_ssl((i_constr-1)*3+3,5)=(coup_nod_ssl(i_constr,2)-1)*3+3; 
    end
end

%rail pad-sleeper
Constr_rps=zeros(n_nod_coupled_rps*dofs,5);% multiply by 3 means 3 dofs
Constr_rps(:,2)=-1;
Constr_rps(:,4)=1;
for i_constr=1:n_nod_coupled_rps
    switch dofs
        case 1
            Constr_rps((i_constr-1)*1+1,3)=coup_nod_rps(i_constr,1);
            Constr_rps((i_constr-1)*1+1,5)=coup_nod_rps(i_constr,2);
        case 2
            Constr_rps((i_constr-1)*2+1,3)=(coup_nod_rps(i_constr,1)-1)*2+1;
            Constr_rps((i_constr-1)*2+1,5)=(coup_nod_rps(i_constr,2)-1)*2+1;
            Constr_rps((i_constr-1)*2+2,3)=(coup_nod_rps(i_constr,1)-1)*2+2;
            Constr_rps((i_constr-1)*2+2,5)=(coup_nod_rps(i_constr,2)-1)*2+2;
        case 3 
            Constr_rps((i_constr-1)*3+1,3)=(coup_nod_rps(i_constr,1)-1)*3+1;
            Constr_rps((i_constr-1)*3+1,5)=(coup_nod_rps(i_constr,2)-1)*3+1;
            Constr_rps((i_constr-1)*3+2,3)=(coup_nod_rps(i_constr,1)-1)*3+2;
            Constr_rps((i_constr-1)*3+2,5)=(coup_nod_rps(i_constr,2)-1)*3+2;
            Constr_rps((i_constr-1)*3+3,3)=(coup_nod_rps(i_constr,1)-1)*3+3;
            Constr_rps((i_constr-1)*3+3,5)=(coup_nod_rps(i_constr,2)-1)*3+3;
    end
end
%rail-rail pad
Constr_rpr=zeros(n_nod_coupled_rpr*dofs,5);% multiply by 3 means 3 dofs
Constr_rpr(:,2)=-1;
Constr_rpr(:,4)=1;
for i_constr=1:n_nod_coupled_rpr
    switch dofs
        case 1
            Constr_rpr((i_constr-1)*1+1,3)=coup_nod_rpr(i_constr,1);
            Constr_rpr((i_constr-1)*1+1,5)=coup_nod_rpr(i_constr,2);
        case 2
            Constr_rpr((i_constr-1)*2+1,3)=(coup_nod_rpr(i_constr,1)-1)*2+1;
            Constr_rpr((i_constr-1)*2+1,5)=(coup_nod_rpr(i_constr,2)-1)*2+1;
            Constr_rpr((i_constr-1)*2+2,3)=(coup_nod_rpr(i_constr,1)-1)*2+2;
            Constr_rpr((i_constr-1)*2+2,5)=(coup_nod_rpr(i_constr,2)-1)*2+2;
        case 3
            Constr_rpr((i_constr-1)*3+1,3)=(coup_nod_rpr(i_constr,1)-1)*3+1;
            Constr_rpr((i_constr-1)*3+1,5)=(coup_nod_rpr(i_constr,2)-1)*3+1;
            Constr_rpr((i_constr-1)*3+2,3)=(coup_nod_rpr(i_constr,1)-1)*3+2;
            Constr_rpr((i_constr-1)*3+2,5)=(coup_nod_rpr(i_constr,2)-1)*3+2;
            Constr_rpr((i_constr-1)*3+3,3)=(coup_nod_rpr(i_constr,1)-1)*3+3;
            Constr_rpr((i_constr-1)*3+3,5)=(coup_nod_rpr(i_constr,2)-1)*3+3;
    end
end

if rail_type==1
%rail-rail
Constr_rr=zeros(n_nod_coupled_rr*dofs,5);
Constr_rr(:,2)=-1;
Constr_rr(:,4)=1;
for i_constr=1:n_nod_coupled_rr
    switch dofs
        case 1
            Constr_rr((i_constr-1)*1+1,3)=coup_nod_rr(i_constr,1);
            Constr_rr((i_constr-1)*1+1,5)=coup_nod_rr(i_constr,2);
        case 2
            Constr_rr((i_constr-1)*2+1,3)=(coup_nod_rr(i_constr,1)-1)*2+1;
            Constr_rr((i_constr-1)*2+1,5)=(coup_nod_rr(i_constr,2)-1)*2+1;
            Constr_rr((i_constr-1)*2+2,3)=(coup_nod_rr(i_constr,1)-1)*2+2;
            Constr_rr((i_constr-1)*2+2,5)=(coup_nod_rr(i_constr,2)-1)*2+2;
        case 3
            Constr_rr((i_constr-1)*3+1,3)=(coup_nod_rr(i_constr,1)-1)*3+1;
            Constr_rr((i_constr-1)*3+1,5)=(coup_nod_rr(i_constr,2)-1)*3+1;
            Constr_rr((i_constr-1)*3+2,3)=(coup_nod_rr(i_constr,1)-1)*3+2;
            Constr_rr((i_constr-1)*3+2,5)=(coup_nod_rr(i_constr,2)-1)*3+2;
            Constr_rr((i_constr-1)*3+3,3)=(coup_nod_rr(i_constr,1)-1)*3+3;
            Constr_rr((i_constr-1)*3+3,5)=(coup_nod_rr(i_constr,2)-1)*3+3;
    end
end
end

Constr=[Constr_ss;Constr_ssl;Constr_rps;Constr_rpr;];%Constr_rr];

figure
plot(Constr(:,3),'ob')
hold on
plot(Constr(:,5),'xr')

T1=zeros(nFeDof);
t1=zeros(nFeDof);
[T1,t1]=addconstrT(Constr,T1,t1);
end