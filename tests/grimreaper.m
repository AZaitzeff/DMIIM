function [u,v,w,vel] = grimreaper(angle1,angle2,T,n,mu)
mu12=(mu(1)+mu(2))/2;
mu13=(mu(1)+mu(3))/2;
angle1r=angle1/180*pi;
angle2r=angle2/180*pi;
angle3r=2*pi-angle1r-angle2r;
d=sin(angle3r)/sin(angle2r)*mu12/mu13;
a=d*(angle2r-pi/2)/(2*(angle3r-pi/2+d*(angle2r-pi/2)));
beta=(pi-angle1r)/(a-a*d+1/2*d);
c=cos(beta*a)^d/cos(beta*d*(1/2-a));

vel=beta*sin(angle3r)*mu12;

u = zeros(n,n);
v = zeros(n,n);
w = ones(n,n);

x = [1:n]/n-1/n-.25;
y = [1:n]/n-1/n;
x2 = x.*(x>=0).*(x<=.5)+(-x).*(x<0)+(1-x).*(x>.5);
y1 = (0.75-1/n) + 1/beta * log(cos(beta*x2)).*(x2<=a)+ 1/(d*beta) * log(c*cos(d*beta*(x2-1/2))).*(x2>a)-vel*T;
y2 = (.25) - (1/beta * log(cos(beta*x2)).*(x2<=a)+ 1/(d*beta) * log(c*cos(d*beta*(x2-1/2))).*(x2>a))+vel*T;
ymax=(0.75-1/n) + 1/beta * log(cos(beta*a));
ymin=(.25) - (1/beta * log(cos(beta*a)));
for i=1:n
  for j=1:n
    yy = y(i);
    xx = x2(j);
    if (yy < y1(j)) && (yy > y2(j))
      u(i,j)= max(yy - y1(j), y2(j)-yy);
      val=min([y1(j)-yy,abs(xx-a),yy-y2(j)]);
      if yy<ymax && yy>ymin
          otherval=-abs(xx-a);
      else
          otherval=-sqrt((xx-a)^2+min((yy-ymax)^2,(yy-ymin)^2));
      end
      if xx<=a
          v(i,j) = val;
          w(i,j) = otherval;
          
      else
          v(i,j) = otherval;
          w(i,j) = val;
          
      end
      
    else
        val=min(abs(yy - y1(j)),abs(yy - y2(j)));
        u(i,j)= val;
        if xx<=a
          v(i,j) = -val;
          w(i,j) = -sqrt((xx-a)^2+min((yy-ymax)^2,(yy-ymin)^2));
          
       else
          v(i,j) = -sqrt((xx-a)^2+min((yy-ymax)^2,(yy-ymin)^2));
          w(i,j) = -val;
          
        end
     end
  end
end