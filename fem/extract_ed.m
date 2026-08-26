  function [U_elem]=extract_ed(edof, U_global)
% ed=extract_ed(edof,a)
%-------------------------------------------------------------
% PURPOSE
%  Extract element displacements from the global displacement
%  vector according to the topology matrix edof.
%
% INPUT:   U:  the global displacement vector
%
%         edof:  topology matrix
%
% OUTPUT: U_elem:  element displacement matrix
%-------------------------------------------------------------

% LAST MODIFIED: H Yingvar 2026-01-13
% Change Log:
% 1. Remove the element number index from the first column of the `edof` matrix.
% 2. Replace loops with matrix-based batch assignment.

%-------------------------------------------------------------
U_elem = U_global(edof')';
%--------------------------end--------------------------------
