function [x1,y1,x2,y2,x3,y3,vel]=VIIM_initialdata_nonsym(angle1,angle2,mu,N,flip,T)

angle1r=angle1/180*pi;
angle2r=angle2/180*pi;
angle3r=2*pi-angle1r-angle2r;
gamma=zeros(1,3);
gamma(1)=(sin(angle3r)+sin(angle2r)-sin(angle1r));
gamma(2)=(sin(angle3r)+sin(angle1r)-sin(angle2r));
gamma(3)=(sin(angle1r)+sin(angle2r)-sin(angle3r));

factor12=(gamma(1)*mu(1)+gamma(2)*mu(2))/2;
factor13=(gamma(1)*mu(1)+gamma(3)*mu(3))/2;
d=factor12/factor13;
a=d*(angle2r-pi/2)/(2*(angle3r-pi/2+d*(angle2r-pi/2)));
beta=(pi-angle1r)/(a-a*d+1/2*d);
c=cos(beta*a)^d/cos(beta*d*(1/2-a));
vel=beta*factor12;

lim1=log(tan(beta*a)+sec(beta*a))/beta;
lim2=log(tan(beta*d*(.5-a))+sec(beta*d*(.5-a)))/(d*beta);
rat1=lim1/(lim1+lim2);
rat2=lim2/(lim1+lim2);
N1=floor(2*N*rat1);
N2=floor(2*N*rat2);
s1=linspace(0,lim1,N1);
s2=linspace(0,lim2,N2);
x1 = 2*atan2(exp(beta*s1)-1,exp(beta*s1)+1)/beta;
y1 = log(cos(beta*x1))*1/(beta)-T*vel;
x2 = 2*atan2(1-exp(d*beta*s2),exp(d*beta*s2)+1)/(d*beta)+1/2;
y2 = log(c*cos(beta*d*(.5-x2)))/(d*beta)-T*vel;

M=floor((flip+y1(end))/sqrt((y1(end)-y1(end-1))^2+(x1(end)-x1(end-1))^2));
x3 = (x1(end)-eps)*ones(1,M);
y3 =linspace(y1(end),-flip,M);
