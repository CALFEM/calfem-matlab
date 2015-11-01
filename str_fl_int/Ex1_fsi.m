% example Ex1_fsi
%----------------------------------------------------------------
% PURPOSE 
%    Fluid-structure interaction, sets up a 2-D fluid model
%
%----------------------------------------------------------------

% LAST MODIFIED: G. Sandberg 1996-06-29
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------


%\label{CB1}
% ----------------------------------------------------------------------
dof_fl=[1:1:189]';
el_no_fl=[1:1:160]';

first_el=[1 2 23 22];  
first_el_row=[];
for i=1:20
  first_el_row=[first_el_row; first_el+(i-1)]; 
end

all_el_fl=[];
for i=1:21:160
  all_el_fl=[all_el_fl; first_el_row+(i-1)]; 
end
edof_fl=[el_no_fl,all_el_fl];

x_coord_first_node_row=[0.0:0.5:10]';
y_coord_first_node_row=zeros(21,1);
coord_first_node_row=[x_coord_first_node_row, y_coord_first_node_row];
coord_fl=[]; 
add=[zeros(21,1), zeros(21,1)+0.5];
for i=1:1:9, 
  coord_fl=[coord_fl; coord_first_node_row+(i-1)*add]; 
end
% ----------------------------------------------------------------------

%\label{CB2}
% ----------------------------------------------------------------------
plotpar=[1,1,0];
[ex_fl,ey_fl]=coordxtr(edof_fl,coord_fl,dof_fl,4);
eldraw2(ex_fl,ey_fl,plotpar,edof_fl(:,1));
% ----------------------------------------------------------------------

%\label{CB3}
% ----------------------------------------------------------------------
c=1500; rho=1000;
ep_fl=[1 c rho 3];

ndof_fl=max(max(dof_fl));
K_fl=zeros(ndof_fl,ndof_fl);
M_fl=K_fl;

for i=1:length(ex_fl)
  [ke_fl,me_fl]=aco2i4d(ex_fl(i,:),ey_fl(i,:),ep_fl);
  K_fl=assem(edof_fl(i,:),K_fl,ke_fl);
  M_fl=assem(edof_fl(i,:),M_fl,me_fl);
end

[La_fl,Egv_fl]=eigen(K_fl,M_fl);
freq_fl=sqrt(La_fl)/2/pi;
% ----------------------------------------------------------------------

%\label{CB4}
% ----------------------------------------------------------------------
clf; modnr=1;
for j=1:3, for i=1:3
 modnr=modnr+1;
 Ed=extract(edof_fl,Egv_fl(:,modnr)); 
 Edabs=abs(Ed); 
 const=max(max(Edabs)); Ed=Ed/const;
 h=fill(ex_fl'+10*1.1*(i-1),ey_fl'-4*1.6*(j-1),Ed');
 MOD= num2str(freq_fl(modnr));
 text(10*1.1*(i-1)+4,-4*1.6*(j-1)-1, MOD)
 set(h,'edgecolor','none')
 hold on, 
end, end
axis equal, axis off
% ----------------------------------------------------------------------

