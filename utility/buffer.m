function [u]=buffer(u,lin)
% BUFFER  pads array with boundary condition.
%   [u]=BUFFER(u,2) mirror boundary condition
%   See also UPBUFFER.
u=padarray(u,[0,2]);
u=upbuffer(u,lin);