% example Ex3_fsi
%----------------------------------------------------------------
% PURPOSE 
%    Fluid-structure interaction, 
%    sets up the complete coupled system 
%    and solves the eigenvalueproblem
%
%----------------------------------------------------------------

% LAST MODIFIED: G. Sandberg 1996-08-07
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

%\label{CB9}
% ----------------------------------------------------------------------
zeroH=0*H;

M=[M_st         zeroH;
   rho*c^2*H'   M_fl];

K=[K_st         -H; 
   zeroH'       K_fl];

[La_p,Egv_p]=eigen(K,M,b);
freq_p=sqrt(real(La_p))/2/pi;
% ----------------------------------------------------------------------


%\label{CB10}
% ----------------------------------------------------------------------
M=[M_st     rho*H;
   zeroH'   M_fl];

K=[K_st     zeroH; 
   -c^2*H'   K_fl];

[La_psi,Egv_psi]=eigen(K,M,b);
freq_psi=sqrt(real(La_psi))/2/pi;
% ----------------------------------------------------------------------


%\label{CB11}
% ----------------------------------------------------------------------
zeroK=0*K_fl;

M=[M_st     zeroH             zeroH;
   zeroH'   rho*c^(-2)*K_fl   zeroK;
   zeroH'   zeroK             zeroK];

K=[K_st     zeroH         -H; 
   zeroH'   zeroK         c^(-2)*K_fl;
   -H'      c^(-2)*K_fl   -(rho*c^2)^(-1)*M_fl];
   
[La_sym,Egv_sym]=eigen(K,M,b);
freq_sym=sqrt(La_sym)/2/pi;
% ----------------------------------------------------------------------


%\label{CB12}
% ----------------------------------------------------------------------
Egv_p_st=Egv_p(1:ndof_st,:);
Egv_p_fl=Egv_p(ndof_st+1:ndof_st+ndof_fl,:);
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



