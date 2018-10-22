function [g,x0,y0]=distancetospline(x,y,x0,pp,pp1,pp2)%unsigned distance
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
g=dist(x,y,x0,y0);
end

function val=dist(x0,y0,x,y)
    val=sqrt((x-x0).^2+(y-y0).^2);
end

function dx=ddist(x0,y0,x,y,y1)
    dx=(x-x0)+(y-y0)*(y1);
end

function ddx=dddist(x0,y0,x,y,y1,y2)
    ddx=1+(y1)^2+(y-y0)*(y2);
end