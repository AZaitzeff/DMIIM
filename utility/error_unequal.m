function [error,trueyl,trueyr]=error_unequal(xl,yl,xr,yr,angle1,angle2,mu,T)

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

trueyl = log(cos(beta*xl))*1/(beta)-T*vel;
trueyr = log(c*cos(beta*d*(.5-xr)))/(d*beta)-T*vel;

nbhrcor=.04;

mask=abs((xl-a))>nbhrcor;
difference=trueyl-yl;
x1=xl(mask);
dif=difference(mask);
total1=(abs(trueyl(1))*(x1(end)));
error1=trap_per(dif,x1);

mask=abs((xr-a))>nbhrcor;
difference=trueyr-yr;
x1=xr(mask);
dif=difference(mask);
error2=trap_per(dif,x1);
total2=(abs(trueyl(1))*(.5-x1(end)));
error=(error1+error2)/(total1+total2);

