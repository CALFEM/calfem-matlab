  function [ed]=extract(edof,a)
% ed=extract(edof,a)
%-------------------------------------------------------------
% PURPOSE
%  Extract element displacements from the global displacement
%  vector according to the topology matrix edof.
%
% INPUT:   a:  the global displacement vector
%
%         edof:  topology matrix
%
% OUTPUT: ed:  element displacement matrix
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa 1993-08-24
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
%-------------------------------------------------------------
    [nie,n]=size(edof);
%
    t=edof(:,2:n);
%
    for i = 1:nie
        ed(i,1:(n-1))=a(t(i,:))';
    end
%
%--------------------------end--------------------------------
