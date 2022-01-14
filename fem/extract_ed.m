  function [ed]=extract_ed(edof,a)
% ed=extract_ed(edof,a)
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
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
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
