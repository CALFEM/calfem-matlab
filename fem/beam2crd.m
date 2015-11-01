    function [excd,eycd]=beam2crd(ex,ey,ed,mag)
%-------------------------------------------------------------
%    PURPOSE
%     Calculate the element continous displacements for a 
%     number of identical 2D Bernoulli beam elements. 
%  
%    INPUT:  ex,ey,
%            ed,
%            mag 
%
%    OUTPUT: excd,eycd 
%-------------------------------------------------------------

% LAST MODIFIED: P-E AUSTRELL 1993-10-15
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%   
  [nie,ned]=size(ed);
  
  for i=1:nie
  
      b=[ex(i,2)-ex(i,1) ey(i,2)-ey(i,1)];   L=sqrt(b*b');   n=b/L;
 
      G=[n(1) n(2)  0    0    0   0;
        -n(2) n(1)  0    0    0   0;
          0    0    1    0    0   0;
          0    0    0   n(1) n(2) 0;
          0    0    0  -n(2) n(1) 0; 
          0    0    0    0    0   1];
   
     d=ed(i,:)';   dl=G*d ;
 %
     xl=[0:L/20:L]';    one=ones(size(xl));
 %
     Cis=[-1 1;
           L 0]/L; ds=[dl(1);dl(4)];
 %          
     ul=([xl one]*Cis*ds)';      
 %          
     Cib=[ 12   6*L    -12  6*L;
          -6*L -4*L^2  6*L -2*L^2;
           0    L^3     0     0;
          L^3     0     0     0]/L^3; 
 %           
     db=[dl(2);dl(3);dl(5);dl(6)];
 %          
     vl=([xl.^3/6 xl.^2/2 xl one]*Cib*db)';
 %   
     cld=[ ul ;
           vl ]; A=[n(1) -n(2);
                    n(2) n(1)];  cd=A*cld;
 %
     xyc=A(:,1)*xl'+ [ex(i,1);ey(i,1)]*one';
 %                 
     excd(i,:)=xyc(1,:)+mag*cd(1,:);   
     eycd(i,:)=xyc(2,:)+mag*cd(2,:);     
  end
%--------------------------end--------------------------------
