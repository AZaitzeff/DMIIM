function [nx,ny]=rotatepts(x,y,theta,x0,y0)
nx=cos(theta)*(x-x0)-sin(theta)*(y-y0)+x0;
ny=sin(theta)*(x-x0)+cos(theta)*(y-y0)+y0;