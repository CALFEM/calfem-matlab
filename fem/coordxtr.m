function [Ex,Ey,Ez]=coordxtr(Edof,Coord,Dof,nen)
%[Ex,Ey,Ez]=coordxtr(Edof,Coord,Dof,nen)
%-------------------------------------------------------------
% PURPOSE
%    Extract nodal coordinate data from the global coordinate 
%    matrix for a number of elements with equal number of 
%    element nodes and dof's. 
%
% INPUT:  Edof :  topology matrix , dim(t)= nie x ned+1
%                         nie= number of identical elements
%                           ned= number of element dof's  
%
%         Coord: global coordinate matrix 
%
%         Dof:   global nodal dof matrix 
%
%         nen:   number of element nodes
%
% OUTPUT: Ex,Ey,Ez : element coordinate matrices
%         Ex=[x1 x2 ...xnen;    one row for each element
%             ...     ...  ;
%             nel     ...  ]  
%             dim= nel x nen ;   nel:number of elemnts 
%-------------------------------------------------------------

% LAST MODIFIED: P-E Austrell 1993-10-14 
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
    [nel,dum]=size(Edof);
     ned=dum-1;
    [n,nsd]=size(Coord);
    [n,nd]=size(Dof);
     nend=ned/nen;
 %
    for i = 1:nel
       nodnum=zeros(1,nen);
       for j = 1:nen
          check=Dof(:,1:nend)-ones(n,1)*Edof(i,(j-1)*nend+2:j*nend+1);
          [indx,dum]=find(check==0);
          nodnum(j)=indx(1);
       end
 %       
          Ex(i,:)=Coord(nodnum,1)';
       if nsd>1
          Ey(i,:)=Coord(nodnum,2)';
       end
       if nsd>2    
          Ez(i,:)=Coord(nodnum,3)';
       end
    end  
%--------------------------end--------------------------------
