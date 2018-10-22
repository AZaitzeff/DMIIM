function [h,hx,hy,hxx,hxy,hyy]=distancetosplineall(x,y,x0,theta,pp,pp1,pp2)
tol=1e-10;
maxiter=5000;
me=1e-12;
for i=1:maxiter
    y0=fnval(pp,x0);
    y01=fnval(pp1,x0);
    y02=fnval(pp2,x0);
    dx=ddist(x,y,x0,y0,y01);
    ddx=dddist(x,y,x0,y0,y01,y02);
    if abs(dx)<tol
        break;
    end
    x0=x0-dx./(ddx+me);
end
h=dist(x,y,x0,y0);
hx=(-cos(theta)*y01+sin(theta))/sqrt(1+y01^2);
hy=(sin(theta)*y01+cos(theta))/sqrt(1+y01^2);
g=1/(1+y01^2+(y0-y)*y02)*(y02/(1+y01^2)^(3/2));
hxx=(-cos(theta)^2-2*sin(theta)*cos(theta)*y01-sin(theta)^2*y01^2)*g;
hxy=(cos(theta)*sin(theta)-cos(2*theta)*y01-sin(theta)*cos(theta)*y01^2)*g;
hyy=(-sin(theta)^2+2*sin(theta)*cos(theta)*y01-cos(theta)^2*y01^2)*g;
end

function val=dist(x0,y0,x,y)
    val=sign(y0-y)*sqrt((x-x0).^2+(y-y0).^2);
end

function dx=ddist(x0,y0,x,y,y1)
    dx=(x-x0)+(y-y0)*(y1);
end

function ddx=dddist(x0,y0,x,y,y1,y2)
    ddx=1+(y1)^2+(y-y0)*(y2);
end