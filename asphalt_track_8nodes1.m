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


function [feNod,feElt,feTyp,feMat,feSec,T1,t1,nodemaster,P]=asphalt_track_8nodes1(width_t,bl_s,h_s,bu_s,z0_s,n1_s,n2_s,h_b,bu_b,bl_b,n2_b,n1_b,n2_a,n1_a,h_a,rho_s,Ezz_s,nu_s,...
    rho_b,Ezz_b,nu_b,rho_a,Ezz_a,nu_a,bl_a,bu_a,rho_usl,Ezz_usl,nu_usl,h_usl,n2_usl,n1_usl,l_usl,rho_sl,Ezz_sl,nu_sl,h_sl,n2_sl,n1_sl,l_sl,...
    Ezz_rp,rho_rp,nu_rp,h_rp,l_rp,E_r,l_r,h_r,rho_r,nu_r,A_r,Ixx_r,Izz_r,dofs,rail_type)


% FINITE ELEMENT MODEL

% subgrade bottom
[feNod_s,feElt_s2]=solidtrack_8nodes_new(bl_s,bu_s,n1_s,h_s,n2_s,width_t,z0_s,l_rp,2);
nNod_s=length(feNod_s);
nElt_s=length(feElt_s2);
feElt_s=ones(nElt_s,12);
feElt_s(:,1)=feElt_s2(1:nElt_s,1);
feElt_s(:,5:12)=feElt_s2(1:nElt_s,2:9);
sevde = feElt_s(:,end)==0;
feElt_s(sevde,2)=2;

plot(feNod_s(:,2),feNod_s(:,3),'or')

% asphalt
z0_a=h_s;
[feNod_a2,feElt_a2]=solidtrack_8nodes_new(bu_a,bl_a,n1_a,h_a,n2_a,width_t,z0_a,l_rp,0);
nNod_a=length(feNod_a2);
nElt_a=length(feElt_a2);
feNod_a=feNod_a2;
feNod_a(:,1)=feNod_a(:,1)+feNod_s(end,1);
feElt_a=ones(nElt_a,12);
feElt_a(:,1)=feElt_a2(1:nElt_a,1)+nElt_s;
feElt_a(:,5:12)=feElt_a2(1:nElt_a,2:9)+nNod_s;
feElt_a(:,2:4)=3;
hold on
plot(feNod_a(:,2),feNod_a(:,3),'xb')
% % ballast
% z0_b=h_s+h_a;
% [feNod_b2,feElt_b2]=solidtrack_8nodes_new(bl_b,bu_b,n1_b,h_b,n2_b,width_t,z0_b,l_rp,0);
% nNod_b=length(feNod_b2);
% nElt_b=length(feElt_b2);
% feNod_b=feNod_b2;
% feNod_b(:,1)=feNod_b(:,1)+nNod_s+nNod_a;
% feElt_b=ones(nElt_b,12);
% feElt_b(:,1)=feElt_b2(1:nElt_b,1)+nElt_s+nElt_a;
% feElt_b(:,5:12)=feElt_b2(1:nElt_b,2:9)+nNod_s+nNod_a;
% feElt_b(:,2:4)=4;
% hold on
% plot(feNod_b(:,2),feNod_b(:,3),'^g')
% %under sleeper pad
% if h_usl==0
%     nNod_usl=0;
%     nElt_usl=0;
%     feNod_usl=[];
%     feElt_usl=[];
% else
%     z0_usl=h_s+h_a+h_b;
%     [feNod_usl2,feElt_usl2]=solidtrack_8nodes_new(l_usl,l_usl,n1_usl,h_usl,n2_usl,width_t,z0_usl,l_rp,0);
%     nNod_usl=length(feNod_usl2);
%     nElt_usl=length(feElt_usl2);
%     feNod_usl=feNod_usl2;
%     feNod_usl(:,1)=feNod_usl(:,1)+nNod_s+nNod_a+nNod_b;
%     feElt_usl=ones(nElt_usl,12);
%     feElt_usl(:,1)=feElt_usl2(1:nElt_usl,1)+nElt_s+nElt_a+nElt_b;
%     feElt_usl(:,5:12)=feElt_usl2(1:nElt_usl,2:9)+nNod_s+nNod_a+nNod_b;
%     feElt_usl(:,2:4)=5;
%     hold on
%     plot(feNod_usl(:,2),feNod_usl(:,3),'om')
% end
% 
% %sleeper
% z0_sl=h_s+h_a+h_b+h_usl;
% [feNod_sl2,feElt_sl2]=solidtrack_8nodes_new(l_sl,l_sl,n1_sl,h_sl,n2_sl,width_t,z0_sl,l_rp,1);
% nNod_sl=length(feNod_sl2);
% nElt_sl=length(feElt_sl2);
% feNod_sl=feNod_sl2;
% feNod_sl(:,1)=feNod_sl(:,1)+nNod_s+nNod_a+nNod_b+nNod_usl;
% feElt_sl=ones(nElt_sl,12);
% feElt_sl(:,1)=feElt_sl2(1:nElt_sl,1)+nElt_s+nElt_a+nElt_b+nElt_usl;
% feElt_sl(:,5:12)=feElt_sl2(1:nElt_sl,2:9)+nNod_s+nNod_a+nNod_b+nNod_usl;
% feElt_sl(:,2:4)=6;
% hold on
% plot(feNod_sl(:,2),feNod_sl(:,3),'om')
% 
% %rail pad
% [~,p]=min(abs(feNod_sl(end-2*n1_sl:end,2)+width_t));
% nNod_rp=16;
% nElt_rp=2;
% feNod_rp(1:nNod_rp,1)=1:1:nNod_rp;
% feNod_rp(1:nNod_rp,1)=feNod_rp(1:nNod_rp,1)+nNod_s+nNod_a+nNod_b+nNod_usl+nNod_sl;
% feNod_rp(1:6,3)=h_s+h_a+h_b+h_sl+h_usl;
% feNod_rp(7:10,3)=h_s+h_a+h_b+h_sl+h_usl+h_rp/2;
% feNod_rp(11:16,3)=h_s+h_a+h_b+h_sl+h_usl+h_rp;
% 
% feNod_rp(1:3,2)=[feNod_sl(end-2*n1_sl+p-1,2)-l_rp/2 feNod_sl(end-2*n1_sl+p-1,2) feNod_sl(end-2*n1_sl+p-1,2)+l_rp/2];
% feNod_rp(4:6,2)=[feNod_sl(end-p+1,2)-l_rp/2 feNod_sl(end-p+1,2) feNod_sl(end-p+1,2)+l_rp/2];
% feNod_rp(7:8,2)=[feNod_sl(end-2*n1_sl+p-1,2)-l_rp/2 feNod_sl(end-2*n1_sl+p-1,2)+l_rp/2];
% feNod_rp(9:10,2)=[feNod_sl(end-p+1,2)-l_rp/2 feNod_sl(end-p+1,2)+l_rp/2];
% feNod_rp(11:13,2)=[feNod_sl(end-2*n1_sl+p-1,2)-l_rp/2 feNod_sl(end-2*n1_sl+p-1,2) feNod_sl(end-2*n1_sl+p-1,2)+l_rp/2];
% feNod_rp(14:16,2)=[feNod_sl(end-p+1,2)-l_rp/2 feNod_sl(end-p+1,2) feNod_sl(end-p+1,2)+l_rp/2];
% 
% feElt_rp=[1+nElt_s+nElt_b+nElt_a+nElt_usl+nElt_sl 7 7 7 feNod_rp(1,1) feNod_rp(3,1) feNod_rp(13,1) feNod_rp(11,1) feNod_rp(2,1) feNod_rp(8,1) feNod_rp(12,1) feNod_rp(7,1);
%     2+nElt_s+nElt_b+nElt_a+nElt_usl+nElt_sl 7 7 7 feNod_rp(4,1) feNod_rp(6,1) feNod_rp(16,1) feNod_rp(14,1) feNod_rp(5,1) feNod_rp(10,1) feNod_rp(15,1) feNod_rp(9,1)];
%     
% hold on
% plot(feNod_rp(:,2),feNod_rp(:,3),'*k')
% 
% %rail
% if rail_type==1
%     feNod_r=feNod_rp;
%     feNod_r(:,1)=feNod_r(:,1)+nNod_rp;
%     feNod_r(1:6,3)=feNod_r(11:16,3);
%     feNod_r(11:16,3)=feNod_r(11:16,3)+h_r;
%     feNod_r(7:10,3)=feNod_r(7:10,3)+h_r/2+h_rp/2;
% 
%     feElt_r=feElt_rp;
%     feElt_r(:,1)=feElt_r(:,1)+nElt_rp;
%     feElt_r(:,2:4)=8;
%     feElt_r(:,5:12)=feElt_r(:,5:12)+nNod_rp;
%     
%     hold on
%     plot(feNod_r(:,2),feNod_r(:,3),'^r')
% elseif rail_type==2
%     feNod_r(1,1)=1+nNod_s+nNod_a+nNod_b+nNod_usl+nNod_sl+nNod_rp;
%     feNod_r(2,1)=2+nNod_s+nNod_a+nNod_b+nNod_usl+nNod_sl+nNod_rp;
%     feNod_r(1,2)=-width_t; %feNod_rp(2,2); 
%     feNod_r(2,2)=+width_t; %feNod_rp(5,2); 
%     feNod_r(:,3)=h_s+h_a+h_b+h_usl+h_sl+h_rp;
%     feElt_r=[1+nElt_s+nElt_b+nElt_a+nElt_usl+nElt_sl+nElt_rp 8 8 8 feNod_r(1,1) 0 0 0 0 0 0 0 ;
%          2+nElt_s+nElt_b+nElt_a+nElt_usl+nElt_sl+nElt_rp 8 8 8 feNod_r(2,1) 0 0 0 0 0 0 0 ];
%     hold on
%     plot(feNod_r(:,2),feNod_r(:,3),'^r')
% end

if rail_type==1 % solid element
    
    feTyp={1 'volu8';2 'volu7';3 'volu8'};
    feSec=[];
    feMat=[1 Ezz_s nu_s rho_s;
    2 Ezz_s nu_s rho_s;
    3 Ezz_a nu_a rho_a];

elseif rail_type==2  % beam element
    
    feTyp={1 'volu7';2 'volu7';3 'volu8';4 'volu8'; 5 'volu8';6 'volu8';7 'volu8';8 'beam1'};
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
        


feElt=[feElt_s; feElt_a];
feNod=[feNod_s; feNod_a];

%LOAD ON RAIL
nFeDof=dofs*size(feNod,1);% depends on dof
P=zeros(nFeDof,1);
if rail_type==1
%     P(end-dofs)=0.5;
%     P(end-3*dofs)=0.5; % depends on the dof, (en3-2)*nof+nof
   P(end-20*dofs)=160e3;

elseif rail_type==2
    P(end)=0.5;
    P(end-dofs)=0.5;
end
% plot(feNod(:,2),feNod(:,3),'or')

%coupled nodes are identified

%soil-track interface
cpind=[1:1:(n1_s+1)]';
indnm=find(feNod(cpind,2)==0);
nodemaster=feNod(indnm,1);
cpind(nodemaster)=[];
n_nod_coupled_ss=n1_s;

%subgrade-asphalt interface
n_nod_coupled_sa=2*n1_a+1;
coup_nod_sa=zeros(n_nod_coupled_sa,2);
coup_nod_sa(:,1)=feNod_a(1:2*n1_a+1,1);
for inode=1:n_nod_coupled_sa
    coup_nod_sa(inode,2)=feNod_s(find((abs(feNod_s(:,2)-feNod(coup_nod_sa(inode,1),2))<1e-10)&(abs(feNod_s(:,3)-feNod(coup_nod_sa(inode,1),3))<1e-10)),1);
end
% %asphalt-ballast interface
% n_nod_coupled_ab=2*n1_b+1;
% coup_nod_ab=zeros(n_nod_coupled_ab,2);
% coup_nod_ab(:,1)=feNod_b(1:2*n1_b+1,1);
% for inode=1:n_nod_coupled_ab
%     coup_nod_ab(inode,2)=feNod_a(find((abs(feNod_a(:,2)-feNod(coup_nod_ab(inode,1),2))<1e-10)&(abs(feNod_a(:,3)-feNod(coup_nod_ab(inode,1),3))<1e-10)),1);
% end
% 
% if h_usl==0
%     %ballast-sleeper pad interface
%     n_nod_coupled_busl=2*n1_usl+1;
%     coup_nod_busl=zeros(n_nod_coupled_busl,2);
%     coup_nod_busl(:,1)=feNod_sl(1:2*n1_usl+1,1);
%     for inode=1:n_nod_coupled_busl
%         coup_nod_busl(inode,2)=feNod_b(find((abs(feNod_b(:,2)-feNod(coup_nod_busl(inode,1),2))<1e-10)&(abs(feNod_b(:,3)-feNod(coup_nod_busl(inode,1),3))<1e-10)),1);
%     end
%     %under sleeper pad-sleeper interface
%     n_nod_coupled_uslsl=0;
%     coup_nod_uslsl=[];
% else
%     %ballast-under sleeper pad interface
%     n_nod_coupled_busl=2*n1_usl+1;
%     coup_nod_busl=zeros(n_nod_coupled_busl,2);
%     coup_nod_busl(:,1)=feNod_usl(1:2*n1_usl+1,1);
%     for inode=1:n_nod_coupled_busl
%         coup_nod_busl(inode,2)=feNod_b(find((abs(feNod_b(:,2)-feNod(coup_nod_busl(inode,1),2))<1e-10)&(abs(feNod_b(:,3)-feNod(coup_nod_busl(inode,1),3))<1e-10)),1);
%     end
%     
%     %under sleeper pad-sleeper interface
%     n_nod_coupled_uslsl=2*n1_sl+1;
%     coup_nod_uslsl=zeros(n_nod_coupled_uslsl,2);
%     coup_nod_uslsl(:,1)=feNod_sl(1:2*n1_sl+1,1);
%     for inode=1:n_nod_coupled_uslsl
%         coup_nod_uslsl(inode,2)=feNod_usl(find((abs(feNod_usl(:,2)-feNod(coup_nod_uslsl(inode,1),2))<1e-10)&(abs(feNod_usl(:,3)-feNod(coup_nod_uslsl(inode,1),3))<1e-10)),1);
%     end
% end
% %rail pad-sleeper
% n_nod_coupled_rps=6;
% coup_nod_rps=zeros(n_nod_coupled_rps,2);
% coup_nod_rps(:,1)=feNod_rp(1:6,1);
%  for inode=1:n_nod_coupled_rps
%         coup_nod_rps(inode,2)=feNod_sl(find((abs(feNod_sl(:,2)-feNod(coup_nod_rps(inode,1),2))<1e-10)&(abs(feNod_sl(:,3)-feNod(coup_nod_rps(inode,1),3))<1e-10)),1);
%  end
% 
% 
% % for inode=1:n_nod_coupled_rps/2
% %     coup_nod_rps(inode,2)=feNod_sl(end-n1_sl+p-1,1);
% % end
% % 
% % for inode=n_nod_coupled_rps/2+1:n_nod_coupled_rps
% %     coup_nod_rps(inode,2)=feNod_sl(end-p+1,1);
% % end
% 
% %rail-rail pad
% if rail_type==1
%     n_nod_coupled_rpr=6;
%     
%     coup_nod_rpr=zeros(n_nod_coupled_rps,2);
%     coup_nod_rpr(:,1)=feNod_r(1:6,1);
%     for inode=1:n_nod_coupled_rpr
%         coup_nod_rpr(inode,2)=feNod_rp(find((abs(feNod_rp(:,2)-feNod(coup_nod_rpr(inode,1),2))<1e-10)&(abs(feNod_rp(:,3)-feNod(coup_nod_rpr(inode,1),3))<1e-10)),1);
%     end
%     
% elseif rail_type==2
%     coup_nod_rpr(1:3,1)=feNod_r(1,1);
%     coup_nod_rpr(1,2)=feNod_rp(11,1);
%     coup_nod_rpr(2,2)=feNod_rp(12,1);
%     coup_nod_rpr(3,2)=feNod_rp(13,1);
%     coup_nod_rpr(4:6,1)=feNod_r(2,1);
%     coup_nod_rpr(4,2)=feNod_rp(14,1);
%     coup_nod_rpr(5,2)=feNod_rp(15,1);
%     coup_nod_rpr(6,2)=feNod_rp(16,1);
%      
%     n_nod_coupled_rpr=6;
% end
% %rail-rail
% if rail_type==1
%     n_nod_coupled_rr=4;
%     coup_nod_rr=zeros(n_nod_coupled_rr,2);
%     nodrail1=[11 13 14 16];
%     nodrail2=[12 12 15 15];
%     coup_nod_rr(:,1)=feNod_r(nodrail1,1);
%     coup_nod_rr(:,2)=feNod_r(nodrail2,1);
% end

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
Constr_sa=zeros(n_nod_coupled_sa*dofs,5); % multiply by 3 means 3 dofs 
Constr_sa(:,2)=-1;
Constr_sa(:,4)=1;
for i_constr=1:n_nod_coupled_sa
    switch dofs
        case 1
            Constr_sa((i_constr-1)*1+1,3)=coup_nod_sa(i_constr,1);
            Constr_sa((i_constr-1)*1+1,5)=coup_nod_sa(i_constr,2);
        case 2
            Constr_sa((i_constr-1)*2+1,3)=(coup_nod_sa(i_constr,1)-1)*2+1;
            Constr_sa((i_constr-1)*2+1,5)=(coup_nod_sa(i_constr,2)-1)*2+1;
            Constr_sa((i_constr-1)*2+2,3)=(coup_nod_sa(i_constr,1)-1)*2+2;
            Constr_sa((i_constr-1)*2+2,5)=(coup_nod_sa(i_constr,2)-1)*2+2;
        case 3           
            Constr_sa((i_constr-1)*3+1,3)=(coup_nod_sa(i_constr,1)-1)*3+1;
            Constr_sa((i_constr-1)*3+1,5)=(coup_nod_sa(i_constr,2)-1)*3+1;% multiply by 3 means 3 dofs
            Constr_sa((i_constr-1)*3+2,3)=(coup_nod_sa(i_constr,1)-1)*3+2;
            Constr_sa((i_constr-1)*3+2,5)=(coup_nod_sa(i_constr,2)-1)*3+2; 
            Constr_sa((i_constr-1)*3+3,3)=(coup_nod_sa(i_constr,1)-1)*3+3;
            Constr_sa((i_constr-1)*3+3,5)=(coup_nod_sa(i_constr,2)-1)*3+3; 
    end
end
% %asphalt-ballast interface
% Constr_ab=zeros(n_nod_coupled_ab*dofs,5);% multiply by 3 means 3 dofs
% Constr_ab(:,2)=-1;
% Constr_ab(:,4)=1;
% for i_constr=1:n_nod_coupled_ab
%     switch dofs
%         case 1
%             Constr_ab((i_constr-1)*1+1,3)=coup_nod_ab(i_constr,1);
%             Constr_ab((i_constr-1)*1+1,5)=coup_nod_ab(i_constr,2);
%         case 2
%             Constr_ab((i_constr-1)*2+1,3)=(coup_nod_ab(i_constr,1)-1)*2+1;
%             Constr_ab((i_constr-1)*2+1,5)=(coup_nod_ab(i_constr,2)-1)*2+1;
%             Constr_ab((i_constr-1)*2+2,3)=(coup_nod_ab(i_constr,1)-1)*2+2;
%             Constr_ab((i_constr-1)*2+2,5)=(coup_nod_ab(i_constr,2)-1)*2+2;
%         case 3    
%             Constr_ab((i_constr-1)*3+1,3)=(coup_nod_ab(i_constr,1)-1)*3+1;
%             Constr_ab((i_constr-1)*3+1,5)=(coup_nod_ab(i_constr,2)-1)*3+1;
%             Constr_ab((i_constr-1)*3+2,3)=(coup_nod_ab(i_constr,1)-1)*3+2;
%             Constr_ab((i_constr-1)*3+2,5)=(coup_nod_ab(i_constr,2)-1)*3+2;
%             Constr_ab((i_constr-1)*3+3,3)=(coup_nod_ab(i_constr,1)-1)*3+3;
%             Constr_ab((i_constr-1)*3+3,5)=(coup_nod_ab(i_constr,2)-1)*3+3;
%     end
% end
% if h_usl==0
%     Constr_busl=zeros(n_nod_coupled_busl*dofs,5);% multiply by 3 means 3 dofs
%     Constr_busl(:,2)=-1;
%     Constr_busl(:,4)=1;
%     for i_constr=1:n_nod_coupled_busl
%         switch dofs
%             case 1
%                 Constr_busl((i_constr-1)*1+1,3)=coup_nod_busl(i_constr,1);
%                 Constr_busl((i_constr-1)*1+1,5)=coup_nod_busl(i_constr,2);
%             case 2
%                 Constr_busl((i_constr-1)*2+1,3)=(coup_nod_busl(i_constr,1)-1)*2+1;
%                 Constr_busl((i_constr-1)*2+1,5)=(coup_nod_busl(i_constr,2)-1)*2+1;
%                 Constr_busl((i_constr-1)*2+2,3)=(coup_nod_busl(i_constr,1)-1)*2+2;
%                 Constr_busl((i_constr-1)*2+2,5)=(coup_nod_busl(i_constr,2)-1)*2+2;
%             case 3                
%                 Constr_busl((i_constr-1)*3+1,3)=(coup_nod_busl(i_constr,1)-1)*3+1;
%                 Constr_busl((i_constr-1)*3+1,5)=(coup_nod_busl(i_constr,2)-1)*3+1;
%                 Constr_busl((i_constr-1)*3+2,3)=(coup_nod_busl(i_constr,1)-1)*3+2;
%                 Constr_busl((i_constr-1)*3+2,5)=(coup_nod_busl(i_constr,2)-1)*3+2;
%                 Constr_busl((i_constr-1)*3+3,3)=(coup_nod_busl(i_constr,1)-1)*3+3;
%                 Constr_busl((i_constr-1)*3+3,5)=(coup_nod_busl(i_constr,2)-1)*3+3;
%         end
%     end
%     Constr_uslsl=[];
% else
%     %ballast-under sleeper pad interface
%     Constr_busl=zeros(n_nod_coupled_busl*3,5);% multiply by 3 means 3 dofs
%     Constr_busl(:,2)=-1;
%     Constr_busl(:,4)=1;
%     for i_constr=1:n_nod_coupled_busl
%         switch dofs
%             case 1
%                 Constr_busl((i_constr-1)*1+1,3)=coup_nod_busl(i_constr,1);
%                 Constr_busl((i_constr-1)*1+1,5)=coup_nod_busl(i_constr,2);
%             case 2
%                 Constr_busl((i_constr-1)*2+1,3)=(coup_nod_busl(i_constr,1)-1)*2+1;
%                 Constr_busl((i_constr-1)*2+1,5)=(coup_nod_busl(i_constr,2)-1)*2+1;
%                 Constr_busl((i_constr-1)*2+2,3)=(coup_nod_busl(i_constr,1)-1)*2+2;
%                 Constr_busl((i_constr-1)*2+2,5)=(coup_nod_busl(i_constr,2)-1)*2+2;
%             case 3                
%                 Constr_busl((i_constr-1)*3+1,3)=(coup_nod_busl(i_constr,1)-1)*3+1;
%                 Constr_busl((i_constr-1)*3+1,5)=(coup_nod_busl(i_constr,2)-1)*3+1;
%                 Constr_busl((i_constr-1)*3+2,3)=(coup_nod_busl(i_constr,1)-1)*3+2;
%                 Constr_busl((i_constr-1)*3+2,5)=(coup_nod_busl(i_constr,2)-1)*3+2;
%                 Constr_busl((i_constr-1)*3+3,3)=(coup_nod_busl(i_constr,1)-1)*3+3;
%                 Constr_busl((i_constr-1)*3+3,5)=(coup_nod_busl(i_constr,2)-1)*3+3;
%         end
%     end
%     %under sleeper pad-sleeper interface
%     Constr_uslsl=zeros(n_nod_coupled_uslsl*dofs,5);% multiply by 3 means 3 dofs
%     Constr_uslsl(:,2)=-1;
%     Constr_uslsl(:,4)=1;
%     for i_constr=1:n_nod_coupled_uslsl
%         switch dofs
%             case 1
%                 Constr_uslsl((i_constr-1)*1+1,3)=coup_nod_uslsl(i_constr,1);
%                 Constr_uslsl((i_constr-1)*1+1,5)=coup_nod_uslsl(i_constr,2);
%             case 2
%                 Constr_uslsl((i_constr-1)*2+1,3)=(coup_nod_uslsl(i_constr,1)-1)*2+1;
%                 Constr_uslsl((i_constr-1)*2+1,5)=(coup_nod_uslsl(i_constr,2)-1)*2+1;
%                 Constr_uslsl((i_constr-1)*2+2,3)=(coup_nod_uslsl(i_constr,1)-1)*2+2;
%                 Constr_uslsl((i_constr-1)*2+2,5)=(coup_nod_uslsl(i_constr,2)-1)*2+2;
%             case 3 
%                 Constr_uslsl((i_constr-1)*3+1,3)=(coup_nod_uslsl(i_constr,1)-1)*3+1;
%                 Constr_uslsl((i_constr-1)*3+1,5)=(coup_nod_uslsl(i_constr,2)-1)*3+1;
%                 Constr_uslsl((i_constr-1)*3+2,3)=(coup_nod_uslsl(i_constr,1)-1)*3+2;
%                 Constr_uslsl((i_constr-1)*3+2,5)=(coup_nod_uslsl(i_constr,2)-1)*3+2;
%                 Constr_uslsl((i_constr-1)*3+3,3)=(coup_nod_uslsl(i_constr,1)-1)*3+3;
%                 Constr_uslsl((i_constr-1)*3+3,5)=(coup_nod_uslsl(i_constr,2)-1)*3+3;
%         end
%     end
% end
% %rail pad-sleeper
% Constr_rps=zeros(n_nod_coupled_rps*dofs,5);% multiply by 3 means 3 dofs
% Constr_rps(:,2)=-1;
% Constr_rps(:,4)=1;
% for i_constr=1:n_nod_coupled_rps
%     switch dofs
%         case 1
%             Constr_rps((i_constr-1)*1+1,3)=coup_nod_rps(i_constr,1);
%             Constr_rps((i_constr-1)*1+1,5)=coup_nod_rps(i_constr,2);
%         case 2
%             Constr_rps((i_constr-1)*2+1,3)=(coup_nod_rps(i_constr,1)-1)*2+1;
%             Constr_rps((i_constr-1)*2+1,5)=(coup_nod_rps(i_constr,2)-1)*2+1;
%             Constr_rps((i_constr-1)*2+2,3)=(coup_nod_rps(i_constr,1)-1)*2+2;
%             Constr_rps((i_constr-1)*2+2,5)=(coup_nod_rps(i_constr,2)-1)*2+2;
%         case 3 
%             Constr_rps((i_constr-1)*3+1,3)=(coup_nod_rps(i_constr,1)-1)*3+1;
%             Constr_rps((i_constr-1)*3+1,5)=(coup_nod_rps(i_constr,2)-1)*3+1;
%             Constr_rps((i_constr-1)*3+2,3)=(coup_nod_rps(i_constr,1)-1)*3+2;
%             Constr_rps((i_constr-1)*3+2,5)=(coup_nod_rps(i_constr,2)-1)*3+2;
%             Constr_rps((i_constr-1)*3+3,3)=(coup_nod_rps(i_constr,1)-1)*3+3;
%             Constr_rps((i_constr-1)*3+3,5)=(coup_nod_rps(i_constr,2)-1)*3+3;
%     end
% end
% %rail-rail pad
% Constr_rpr=zeros(n_nod_coupled_rpr*dofs,5);% multiply by 3 means 3 dofs
% Constr_rpr(:,2)=-1;
% Constr_rpr(:,4)=1;
% for i_constr=1:n_nod_coupled_rpr
%     switch dofs
%         case 1
%             Constr_rpr((i_constr-1)*1+1,3)=coup_nod_rpr(i_constr,1);
%             Constr_rpr((i_constr-1)*1+1,5)=coup_nod_rpr(i_constr,2);
%         case 2
%             Constr_rpr((i_constr-1)*2+1,3)=(coup_nod_rpr(i_constr,1)-1)*2+1;
%             Constr_rpr((i_constr-1)*2+1,5)=(coup_nod_rpr(i_constr,2)-1)*2+1;
%             Constr_rpr((i_constr-1)*2+2,3)=(coup_nod_rpr(i_constr,1)-1)*2+2;
%             Constr_rpr((i_constr-1)*2+2,5)=(coup_nod_rpr(i_constr,2)-1)*2+2;
%         case 3
%             Constr_rpr((i_constr-1)*3+1,3)=(coup_nod_rpr(i_constr,1)-1)*3+1;
%             Constr_rpr((i_constr-1)*3+1,5)=(coup_nod_rpr(i_constr,2)-1)*3+1;
%             Constr_rpr((i_constr-1)*3+2,3)=(coup_nod_rpr(i_constr,1)-1)*3+2;
%             Constr_rpr((i_constr-1)*3+2,5)=(coup_nod_rpr(i_constr,2)-1)*3+2;
%             Constr_rpr((i_constr-1)*3+3,3)=(coup_nod_rpr(i_constr,1)-1)*3+3;
%             Constr_rpr((i_constr-1)*3+3,5)=(coup_nod_rpr(i_constr,2)-1)*3+3;
%     end
% end
% 
% if rail_type==1
% %rail-rail
% Constr_rr=zeros(n_nod_coupled_rr*dofs,5);
% Constr_rr(:,2)=-1;
% Constr_rr(:,4)=1;
% for i_constr=1:n_nod_coupled_rr
%     switch dofs
%         case 1
%             Constr_rr((i_constr-1)*1+1,3)=coup_nod_rr(i_constr,1);
%             Constr_rr((i_constr-1)*1+1,5)=coup_nod_rr(i_constr,2);
%         case 2
%             Constr_rr((i_constr-1)*2+1,3)=(coup_nod_rr(i_constr,1)-1)*2+1;
%             Constr_rr((i_constr-1)*2+1,5)=(coup_nod_rr(i_constr,2)-1)*2+1;
%             Constr_rr((i_constr-1)*2+2,3)=(coup_nod_rr(i_constr,1)-1)*2+2;
%             Constr_rr((i_constr-1)*2+2,5)=(coup_nod_rr(i_constr,2)-1)*2+2;
%         case 3
%             Constr_rr((i_constr-1)*3+1,3)=(coup_nod_rr(i_constr,1)-1)*3+1;
%             Constr_rr((i_constr-1)*3+1,5)=(coup_nod_rr(i_constr,2)-1)*3+1;
%             Constr_rr((i_constr-1)*3+2,3)=(coup_nod_rr(i_constr,1)-1)*3+2;
%             Constr_rr((i_constr-1)*3+2,5)=(coup_nod_rr(i_constr,2)-1)*3+2;
%             Constr_rr((i_constr-1)*3+3,3)=(coup_nod_rr(i_constr,1)-1)*3+3;
%             Constr_rr((i_constr-1)*3+3,5)=(coup_nod_rr(i_constr,2)-1)*3+3;
%     end
% end
% end

% Constr=[Constr_ss;Constr_sa;Constr_ab;Constr_busl;Constr_uslsl;Constr_rps;Constr_rpr;];%Constr_rr];
Constr=[Constr_ss; Constr_sa];
% figure
% plot(Constr(:,3),'ob')
% hold on
% plot(Constr(:,5),'xr')

T1=zeros(nFeDof);
t1=zeros(nFeDof);
[T1,t1]=addconstrT(Constr,T1,t1);
end