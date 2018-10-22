function [gx,gy,gx0x0,gx0x,gx0y,x0x,x0y]=distderpartsign(y01,y02,ddx)

normalnorm=sqrt(1+y01^2);

gx=y01/normalnorm;
gy=-1/normalnorm;
gx0x0=-(ddx)*y02/normalnorm^3;
gx0x=y02/normalnorm^3;
gx0y=y01*y02/normalnorm^3;
x0x=1/(ddx);
x0y=y01/(ddx);