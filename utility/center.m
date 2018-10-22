function [du]=center(u)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Backward difference operator 

% matrix u, which is assumed to be buffered (layer thickness 2).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uses: upbuffer2.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[~,n]=size(u);  % Buffered size.

n=n-4;      % Find the unbuffered size.

     % Find the unbuffered size.
du=u;       % Initialization for dx.

% Carry out operation on unbuffered part of u:

du(2:n+3)=( u(3:n+4) - u(1:n+2) )/2;
%du=upbuffer(du);