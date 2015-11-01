function arrow2(x0,y0,la,fi,loc,s)
%-------------------------------------------------------------
% PURPOSE
%  Draw a 2D arrow with length la and angle fi 
%  at location (x0,y0) 
%
% INPUT:  
%  loc=-1,0 or 1 depending on if the tail, the midpoint or 
%                  the tip should be in point (x0,y0) 
%                                      
% OUTPUT: 
%-------------------------------------------------------------

%    REFERENCES
%        P-E AUSTRELL 1994-01-06 
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
%-------------------------------------------------------------
  xyl0=[-0.5 0.35 0.35 0.5 0.35 0.35;
         0   0  -0.05  0  0.05  0  ];
  nar=length(x0);       

  for i=1:nar
     xyl=xyl0-0.5*loc(i)*[1 1 1 1 1 1;
                          0 0 0 0 0 0];
     v=fi(i); L=la(i);

     xy=L*[ cos(v)  -sin(v) ;
            sin(v)   cos(v) ]*xyl+[x0(i)*ones(1,6);
                                   y0(i)*ones(1,6)];
     hold on
     plot( xy(1,:)',xy(2,:)',s)
     hold off
  end
%--------------------------end--------------------------------
