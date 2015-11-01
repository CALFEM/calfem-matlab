 function [H]=assem_ns(edof_s,edof_f,H,He)
% [H]=assem_ns(edof_s,edof_f,H,He)
%                Lund Institute of Technology
%-------------------------------------------------------------
    [nies,ns]=size(edof_s)  
    [nief,nf]=size(edof_f);
    ts=edof_s(:,2:ns)
    tf=edof_f(:,2:nf);
    for i = 1:nies
		H(ts(i,:),tf(i,:)) = H(ts(i,:),tf(i,:))+He;
    end
%--------------------------end--------------------------------
