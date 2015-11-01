  function [Resp,Freq]=freqresp(d,dt)
% [Resp,Freq]=freqresp(d,dt)
%----------------------------------------------------------
% PURPOSE
%  Calculate the fourier transform of a function
%  and plots the result in a loglog scale.
%
% INPUT:
%    d : time history function to be transformed
%    dt: time step used in the creation of d
%
% OUTPUT:
%    Resp: response in frequncy domain
%    Freq: corresponding frequencies
%
%----------------------------------------------------------

% LAST MODIFIED: G Sandberg  1993-11-07
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%----------------------------------------------------------
nr=length(d);
if (nr<=80);             f2=64;  end
if (nr>80  & nr<=150);   f2=128;  end
if (nr>150 & nr<=350);   f2=256;  end
if (nr>350 & nr<=850);   f2=512;  end 
if (nr>850 & nr<=1550);  f2=1024; end 
if (nr>1550);            f2=2048; end 
ff2=f2/2;
Y=fft(d,f2);
Py=Y.*conj(Y)/f2;
Resp=Py(1:ff2);
Freq=(1/dt)*(0:ff2-1)/f2;

%--------------------------end--------------------------------
