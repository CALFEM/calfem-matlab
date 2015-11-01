function Functest(functype)
%function Functest(functype)
%######################## Functest.m ########################  
%------------------------------------------------------------
% PURPOSE
%   Testing of the functions in Calfem.
%   The functions are divided into the follwing groups:
%
%   * static system functions             :functype='statsys'
%   * dynamic system and other functions  :functype='dynsys '
%   * bar element functions               :functype='barel  '
%   * beam element functions              :functype='beamel '
%   * field element functions             :functype='fieldel'
%   * solid element functions             :functype='solidel'
%   * plate bending element functions     :functype='platel '
%   * error handling and other functions  :functype='errorh '

% LAST MODIFIED: Kent Persson 1997-04-15
%------------------------------------------------------------
% 
%           Function type is assigned here:
%
% functype='statsys' 
% functype='dynsys '
% functype='barel  '
% functype='beamel '
% functype='fieldel'
% functype='solidel'
% functype='platel '
% functype='errorh '
% functype='testall' 
%
%       static system functions  :
%
 if (functype=='statsys') | (functype=='testall');
    disp('**************************************')
    disp('***** testing of system functions ****')
    disp('****************************** *******')
    En1=[1 1 2]; 
    En2=[2 2 3];
    K=zeros(3,3);
    f=zeros(3,1); 
    Ke=[1 -1;-1 1];
    fe=[1;1];
    En1=[1 1 2]; En2=[2 2 3];
    disp('[K,f]=assem(En1,K,Ke,f,fe)')
    [K,f]=assem(En1,K,Ke,f,fe)
    [K,f]=assem(En2,K,Ke,f,fe)
    pause
    Bc=[1 0];
    disp('[d,Q]=solveq(K,f,Bc)')
    [d,Q]=solveq(K,f,Bc)
    pause
    disp('ed1=extract(En1,d)')
    ed1=extract(En1,d)
    M=ones(3,3)+eye(3,3)
    b=1;
    disp('[L,X]=eigen(K,M,b)')
    [L,X]=eigen(K,M,b)
    pause
    Cd=1;
    disp('[K1,f1]=statcon(K,f,Cd)')
    [K1,f1]=statcon(K,f,Cd)
    pause
    disp('K1=red(K,Cd)')
    K1=red(K,Cd)
  end;
%
%       dynamic system and other functions  :
%
 if (functype=='dynsys ') | (functype=='testall');
    disp('******************************************')
    disp('** testing of dynamic system functions ***')
    disp('******************************************')
 end;
%................... cont........................... 
    
%
%       bar element functions  :
%
 if (functype=='barel  ') | (functype=='testall');
    disp('**************************************')
    disp('****** testing of bar elements *******')
    disp('****************************** *******') 
    ep=[1];
    disp('Ke=bar1e(ep)')
    Ke=bar1e(ep)
    pause
    ed=[1 3];
    disp('P=bar1s(ep,ed)')
    P=bar1s(ep,ed)
    pause
    ex=[1 2];
    ey=[1 2];
    ep=[1 1];
    disp('Ke=bar2e(ex,ey,ep)')
    Ke=bar2e(ex,ey,ep)
    pause
    ed=[1 1 3 3];
    disp('P=bar2s(ex,ey,ep,ed)')
    P=bar2s(ex,ey,ep,ed)
    pause
    N=1;
    disp('Ke=bar2g(ex,ey,ep,N)')
    Ke=bar2g(ex,ey,ep,N)
    pause
    disp('P=bar2s(ex,ey,ep,ed)')
    P=bar2s(ex,ey,ep,ed)
    pause
    ex=[1 2];
    ey=[1 2];
    ez=[1 2];
    ed=[1 1 1 3 3 3];
    disp('Ke=bar3e(ex,ey,ez,ep)')
    Ke=bar3e(ex,ey,ez,ep)
    pause
    disp('P=bar3s(ex,ey,ez,ep,ed)')
    P=bar3s(ex,ey,ez,ep,ed)
 end;
%
%       beam element functions  :
%
 if (functype=='beamel ')  | (functype=='testall');
    disp('**************************************')
    disp('****** testing of beam elements  *****')
    disp('****************************** *******') 
    ex=[1 2];
    ey=[1 2];
    ep=[ 1 1 1 1 1 1 ];
    disp('[Ke,Me,Ce]=beam2d(ex,ey,ep)')
    [Ke,Me,Ce]=beam2d(ex,ey,ep)
    pause
%-------------------------------------------------------
    ed=[1 1 1 3 3 3]; 
    ev=[1 1 1 3 3 3]; 
    ea=[1 1 1 3 3 3];
    disp('P=beam2ds(ex,ey,ep,ed,ev,ea)')
    P=beam2ds(ex,ey,ep,ed,ev,ea)
    pause
%-------------------------------------------------------
    ep=[1 1 1];
    eq=[1 1];
    disp('[Ke,fe]=beam2e(ex,ey,ep,eq)')
    [Ke,fe]=beam2e(ex,ey,ep,eq)
    pause
%-------------------------------------------------------
    disp('P=beam2s(ex,ey,ep,ed,eq,10)')
    P=beam2s(ex,ey,ep,ed,eq,10)
    pause
%-------------------------------------------------------    
    N=-1;
    q=1;
    disp('[Ke,fe]=beam2g(ex,ey,ep,N,q)')   
    [Ke,fe]=beam2g(ex,ey,ep,N,q)
    pause
%-------------------------------------------------------
    disp('P=beam2gs(ex,ey,ep,ed,N,q)')
    P=beam2gs(ex,ey,ep,ed,N,q)
    pause
%-------------------------------------------------------
    ep=[1 1 1 1 1];
    disp('[Ke,fe]=beam2w(ex,ey,ep,q)')
    [Ke,fe]=beam2w(ex,ey,ep,eq)
    pause
%-------------------------------------------------------
    ex=[0 2];
    ey=[0 2];
    ez=[0 2];
    e0=[-1 0 1];
    ep=[1 1 1 1 1 1];
    disp('Ke=beam3e(ex,ey,ez,e0,ep)')
    Ke=beam3e(ex,ey,ez,e0,ep)
    pause
%-------------------------------------------------------
    ed=[1 1 1 1 1 1 3 3 3 3 3 3];
    disp('P=beam3s(ex,ey,ez,e0,ep,ed)')
    P=beam3s(ex,ey,ez,e0,ep,ed)
  end;
%
%      field element  functions  :
% 
 if (functype=='fieldel')  | (functype=='testall');
    disp('***************************************')
    disp('****** testing of field elements  *****')
    disp('***************************************')
    D=[1 0;
       0 1];
    Q=1;
    ep=2;
    ex=[1 2 2];
    ey=[1 1 2];
    ed=[0 1 1];
    disp('[Ke,fe]=flw2te(ex,ey,ep,D,Q)')
    %Ke=flw2te(ex,ey,D)
    [Ke,fe]=flw2te(ex,ey,ep,D,Q)
    pause
%-------------------------------------------------------
    disp('[es,et]=flw2ts(ex,ey,D,ed)')
    [es,et]=flw2ts(ex,ey,D,ed)
    pause
%-------------------------------------------------------
%    ep=2;
%    disp('[Ke,Ce,fe]=flw2td([ex;ey],ep,D,Q)')
%    [Ke,Ce,fe]=flw2td([ex;ey],ep,D,Q)
%    pause
%-------------------------------------------------------
    ex=[1 2 2 1];
    ey=[1 1 2 2];
    ep=[2,2];
    ed=[0 1 1 0];
    disp('[Ke,fe]=flw2i4e(ex,ey,ep,D,Q)')
    [Ke,fe]=flw2i4e(ex,ey,ep,D,Q)
    pause
%-------------------------------------------------------
    disp('[es,X]=flw2i4s(ex,ey,ep,D,ed)')
    [es,et,xi]=flw2i4s(ex,ey,ep,D,ed) 
    pause
%-------------------------------------------------------
    disp('[Ke,fe]=flw2qe(ex,ey,[1],D,Q)')
    [Ke,fe]=flw2qe(ex,ey,[1],D,Q)
    pause
%-------------------------------------------------------
    disp('[es,et]=flw2qs(ex,ey,ep,D,ed,Q)')
    [es,et]=flw2qs(ex,ey,[1],D,ed,Q)
    pause
%-------------------------------------------------------
    ex=[1 2 2 1 1.5 2 1.5 1 ];
    ey=[1 1 2 2  1 1.5 2 1.5];
    ed=[0 1 1 0 0.5 1 0.5 0 ];
    disp('[Ke,fe]=flw2i8e(ex,ey,ep,D,Q)')
    [Ke,fe]=flw2i8e(ex,ey,ep,D,Q)
    pause
%-------------------------------------------------------
    disp('[es,X]=flw2i8s(ex,ey,ep,D,ed)')
    [es,et,xi]=flw2i8s(ex,ey,ep,D,ed) 
    pause
%------------------------------------------------------- 
    ep=[2];
    Q=[1];
    ex=[1 2 2 1 1 2 2 1];
    ey=[1 1 2 2 1 1 2 2];
    ez=[1 1 1 1 2 2 2 2];
    D=[1 0 0;
       0 1 0; 
       0 0 1];
    ed=[0 1 1 0 0 1 1 0 ];
    disp('[Ke,Fe]=flw3i8e(ex,ey,ez,ep,D,Q)')
    [Ke,Fe]=flw3i8e(ex,ey,ez,ep,D,Q)
    pause
%-------------------------------------------------------
    disp('[es,et,xi]=flw3i8s(ex,ey,ez,ep,D,ed)')
    [es,et,xi]=flw3i8s(ex,ey,ez,ep,D,ed)
  end;
%
%      solid element functions  :
%    
 if (functype=='solidel')  | (functype=='testall');
    disp('***************************************')
    disp('****** testing of solid elements  *****')
    disp('***************************************')
    D=hooke(1,1,0.3);
    b=[1 1]';
    ex=[1 2 2];
    ey=[1 1 2];
    ep=[1 1];
    ed=[0 0 1 0 1 0];
    disp('[Ke,fe]=plante(ex,ey,ep,D,b)')    
    [Ke,fe]=plante(ex,ey,ep,D,b)
    pause
%-------------------------------------------------------
    disp('es=plants(ex,ey,ep,D,ed)')   
    es=plants(ex,ey,ep,D,ed )
    pause
%-------------------------------------------------------
    ex=[1 2 2 1];
    ey=[1 1 2 2];
    ed=[0 0 1 0 1 0 0 0];
    disp('[Ke,fe]=planqe(ex,ey,ep,D,b)')
    [Ke,fe]=planqe(ex,ey,ep,D,b)
    pause
%-------------------------------------------------------
    disp('es=planqs(ex,ey,ep,D,ed,b)')
    es=planqs(ex,ey,ep,D,ed,b)
    pause
%-------------------------------------------------------
    ex=[1 2 ];
    ey=[1 2 ];
    disp('[Ke,fe]=planre(ex,ey,ep,D,b)')
    [Ke,fe]=planre(ex,ey,ep,D,b)
    pause
%-------------------------------------------------------
    disp('es=planrs(ex,ey,ep,D,ed)')
    es=planrs(ex,ey,ep,D,ed)
    pause
%-------------------------------------------------------
    ep=[1 1 1 0.3];
    disp('[Ke,fe]=plantce(ex,ey,ep,b)')
    [Ke,fe]=plantce(ex,ey,ep,b)
    pause
%-------------------------------------------------------
    disp('es=plantcs(ex,ey,ep,ed)')
    es=plantcs(ex,ey,ep,ed)
    pause
%-------------------------------------------------------    
    ex=[1 2 2 1];
    ey=[1 1 2 2];
    ep=[1 1 2];
    disp('[Ke,Fe]=plani4e(ex,ey,ep,D,b)')
    [Ke,Fe]=plani4e(ex,ey,ep,D,b)
    pause
%-------------------------------------------------------
    disp('[es,X]=plani4s(ex,ey,ep,D,ed)')
    [es,X]=plani4s(ex,ey,ep,D,ed)
    pause
%-------------------------------------------------------
    ex=[1 2 2 1 1.5 2 1.5 1 ];
    ey=[1 1 2 2  1 1.5 2 1.5];
    ed=[0 0 1 0 1 0 0 0 0.5 0 1 0 0.5 0 0 0];
    disp('[Ke,Fe]=plani8e(ex,ey,ep,D,b)')
    [Ke,Fe]=plani8e(ex,ey,ep,D,b)
    pause
%-------------------------------------------------------
    disp('[es,X]=plani8s(ex,ey,ep,D,ed)')
    [es,X]=plani8s(ex,ey,ep,D,ed)
    pause
%-------------------------------------------------------
    ex=[1 2 2 1 1 2 2 1];
    ey=[1 1 2 2 1 1 2 2];
    ez=[1 1 1 1 2 2 2 2];
    ep=2;    
    D=hooke(4,1,0.3);
    b=[1 1 1]';
    ed=[0 0 1 0 1 0 0 0 0.5 0 1 0 0.5 0 0 0 0 0 0 0 0 0 0 0];
    disp('[Ke,fe]=soli8e(ex,ey,ez,ep,D,b)')
    [Ke,fe]=soli8e(ex,ey,ez,ep,D,b)
    pause
%-------------------------------------------------------
    disp('[es,X]=soli8s(ex,ey,ez,ep,D,ed)')
    [es,X]=soli8s(ex,ey,ez,ep,D,ed)
%-------------------------------------------------------
  end;
%
%      plate element  functions  :
%  
 if (functype=='platel ')  | (functype=='testall');
    disp('***************************************')
    disp('****** testing of plate element   *****')
    disp('***************************************')
    disp('     TEST FUNCTION NOT IMPLEMENTED     ');
 end;
%
%      error handling and other functions  :
%  
 if (functype=='errorh ')  | (functype=='testall');
 
 end;
