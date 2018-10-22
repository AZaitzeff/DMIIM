function [x1,y1,x2,y2,x3,y3]=VIIM_initialdata_general(angle1,angle2,N,mu12,mu13,mu23,flip,T)
%mu12=(mu(1)+mu(2))/2;
%mu13=(mu(1)+mu(3))/2;
angle1r=angle1/180*pi;
angle2r=angle2/180*pi;
angle3r=2*pi-angle1r-angle2r;
d=sin(angle3r)/sin(angle2r)*mu12/mu13;
a=d*(angle2r-pi/2)/(2*(angle3r-pi/2+d*(angle2r-pi/2)));
beta=(pi-angle1r)/(a-a*d+1/2*d);
c=cos(beta*a)^d/cos(beta*d*(1/2-a));

vel=beta*sin(angle3r)*mu12;

lim1=log(tan(a*beta)+sec(beta*a))/beta;
lim2=log(tan(beta*d*(.5-a))+sec(beta*d*(.5-a)))/(beta*d);
rat1=lim1/(lim1+lim2);
rat2=lim2/(lim1+lim2);
N1=floor(2*N*rat1);
N2=floor(2*N*rat2);
s1=linspace(0,lim1,N1);
s2=linspace(0,lim2,N2);
x1 = 2*atan2(exp(beta*s1)-1,exp(beta*s1)+1)/beta;
y1 = log(cos(beta*x1))*1/(beta)-T*vel;
x2 = 2*atan2(1-exp(beta*d*s2),exp(beta*d*s2)+1)/(beta*d)+1/2;
y2 = log(c*cos(beta*d*(.5-x2)))*1/(d*beta)-T*vel;

M=floor((flip+y1(end))/sqrt((y1(end)-y1(end-1))^2+(x1(end)-x1(end-1))^2));
x3 = (x1(end)-eps)*ones(1,M);
y3 =linspace(y1(end),-flip,M);
