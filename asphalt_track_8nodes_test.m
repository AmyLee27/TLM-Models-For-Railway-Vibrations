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


function [feNod,feElt,feTyp,feMat,feSec,P]=asphalt_track_8nodes_test(width_t,bl_s,h_s,bu_s,z0_s,n1_s,n2_s,h_b,bu_b,bl_b,n2_b,n1_b,n2_a,n1_a,h_a,rho_s,Ezz_s,nu_s,...
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

feElt=[feElt_s];
feNod=[feNod_s];

if rail_type==1 % solid element
    
    feTyp={1 'volu8'};
    feSec=[];
    %feMat=[1 Ezz_s nu_s rho_s];
    feMat=trackMat(feElt,Ezz_s,nu_s,rho_s,Ezz_a,nu_a,rho_a,Ezz_b,nu_b,rho_b,Ezz_usl,nu_usl,rho_usl,Ezz_sl,nu_sl,rho_sl,Ezz_rp,nu_rp,rho_rp,E_r,nu_r,rho_r);

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
        

%LOAD ON RAIL
nFeDof=dofs*size(feNod,1);% depends on dof
P=zeros(nFeDof,1);
if rail_type==1
%     P(end-10*dofs)=1;
%     P(end-30*dofs)=1; % depends on the dof, (en3-2)*nof+nof

    P(end-20*dofs) = 160e3;
    
elseif rail_type==2
    P(end)=0.5;
    P(end-dofs)=0.5;
end
% plot(feNod(:,2),feNod(:,3),'or')

end