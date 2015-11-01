% example Ex2_fsi
%----------------------------------------------------------------
% PURPOSE 
%    Fluid-structure interaction, sets up the structural part
%    and the coupling part.
%
%----------------------------------------------------------------

% LAST MODIFIED: G. Sandberg 1996-08-07
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------


%\label{CB5}
% ----------------------------------------------------------------------
el_no_st=[1:1:20]';
coord_st=[[0:0.5:10]',4*ones(21,1)];

first_dof_st=[1:3:61]';
dof_st=[first_dof_st, first_dof_st+1, first_dof_st+2];

first_edof_st=[1 2 3 4 5 6];
all_el_st=[];
for i=1:1:20, all_el_st=[all_el_st; first_edof_st+3*(i-1)]; end
edof_st=[el_no_st, all_el_st];

plotpar=[1,1,0];
[ex_st,ey_st]=coordxtr(edof_st,coord_st,dof_st,2);
eldraw2(ex_st,ey_st,plotpar,edof_st(:,1));
% ----------------------------------------------------------------------

%\label{CB6}
% ----------------------------------------------------------------------
ndof_st=max(max(dof_st));
E=2.1e11;  A=0.02;  I=1.59e-4;  m=A*2500;
ep_st=[E A I m];
K_st=zeros(ndof_st,ndof_st);       M_st=K_st;

for i=1:length(ex_st)
  [ke_st,me_st]=beam2d(ex_st(i,:),ey_st(i,:),ep_st);
  K_st=assem(edof_st(i,:),K_st,ke_st);
  M_st=assem(edof_st(i,:),M_st,me_st);
end

b=[1 2 61 62];
[La_st,Egv_st]=eigen(K_st,M_st,b);
freq_st=sqrt(La_st)/2/pi;
% ----------------------------------------------------------------------

%\label{CB7}
% ----------------------------------------------------------------------
edof_coup_st=edof_st;
el_no_coup_fl=[141:1:160]';

first_edof_coup_fl=[169 170];
all_el_coup_fl=[];
for i=1:1:20
  all_el_coup_fl=[all_el_coup_fl; first_edof_coup_fl+(i-1)]; 
end
edof_coup_fl=[el_no_coup_fl, all_el_coup_fl];
% ----------------------------------------------------------------------

%%\label{CB8}
% ----------------------------------------------------------------------
H=zeros(ndof_st,ndof_fl);

for i=1:length(ex_st)
  he=cp2s2f(ex_st(i,:),ey_st(i,:),[1]);
  H=assem_ns(edof_coup_st(i,:),edof_coup_fl(i,:),H,he);
end
% ----------------------------------------------------------------------











