% example Ex3_fsi
%----------------------------------------------------------------
% PURPOSE 
%    Fluid-structure interaction, 
%    sets up the modal formulation and solve the
%    eigenvalueproblem
%
%----------------------------------------------------------------

% LAST MODIFIED: G. Sandberg 1996-08-07
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

%\label{CB13}
% ----------------------------------------------------------------------
fp=[rho c]
list_st=[1:1:15];
list_fl=[1:1:50];

[Amodc_pr] = ... 
 fsi_mod(La_st,Egv_st,list_st,La_fl,Egv_fl,list_fl,H,fp,'right');
[La_pr_mod,Egv_pr_mod] = ... 
 eigen(Amodc_pr,eye(length(list_fl)+length(list_st)));
freq_pr_mod=sqrt(real(La_pr_mod))/2/pi;
% ----------------------------------------------------------------------


%\label{CB14}
% ----------------------------------------------------------------------
[Amodc_pl] = ... 
 fsi_mod(La_st,Egv_st,list_st,La_fl,Egv_fl,list_fl,H,fp,'left');
[La_pl_mod,Egv_pl_mod] = ... 
 eigen(Amodc_pl,eye(length(list_fl)+length(list_st)));
freq_pl_mod=sqrt(real(La_pl_mod))/2/pi;
% ----------------------------------------------------------------------


%\label{CB15}
% ----------------------------------------------------------------------
[Egv_pr_red] = ... 
 fsi_egv(Egv_pr_mod,La_st,Egv_st,list_st,La_fl,Egv_fl,list_fl,fp,'right','pr');
% ----------------------------------------------------------------------


%\label{CB16}
% ----------------------------------------------------------------------
Egv_p_st=Egv_p_red(1:ndof_st,:);
Egv_p_fl=Egv_p_red(ndof_st+1:ndof_st+ndof_fl,:);
clf; modnr=1;  
for j=1:3, for i=1:3
 hold on;  
 modnr=modnr+1;
 Ed_st=real(extract(edof_st,Egv_p_st(:,modnr))); 
 Edabs_st=abs(Ed_st); 
 const_st=max(max(Edabs_st)); Ed_st=Ed_st/const_st;
 eldisp2(ex_st+10*1.1*(i-1),ey_st-4*2.5*(j-1)+2.0,Ed_st,[1,1,0],2);
 eldraw2(ex_st+10*1.1*(i-1),ey_st-4*2.5*(j-1)+2.0,[3,1,0]);
 hold on;  
 Ed_fl=real(extract(edof_fl,Egv_p_fl(:,modnr))); 
 Edabs_fl=abs(Ed_fl); 
 const_fl=max(max(Edabs_fl)); Ed_fl=Ed_fl/const_fl;
 h=fill(ex_fl'+10*1.1*(i-1),ey_fl'-4*2.5*(j-1),Ed_fl');
 MOD= num2str(freq_p(modnr));
 text(10*1.1*(i-1)+4,-4*2.5*(j-1)-1, MOD)
 set(h,'edgecolor','none')
end, end
axis equal, axis off
% ----------------------------------------------------------------------




\end{appendix}
