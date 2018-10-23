function [nx,ny,FLAG]=voronoi_x_spline(x1,y1,x3,y3,nx,yrotate,angle,ST)
FLAG=0;
n=numel(nx);

bx1=x1(1:n);
by1=y1(1:n);
[bx1,by1]=buffer1d(bx1,by1,n/8,0);
ppb1 = csapi(bx1,by1);
ppb1d=fnder(ppb1,1);
ppb1dd=fnder(ppb1,2);
%theta=2*pi/3;
theta=(angle/2+(360-angle)/4)/180*pi;
xrotate=.25;
[rx1,ry1]=rotatepts(x1,y1,-theta,xrotate,yrotate);
ppr1 = csapi(rx1,ry1);
ppr1d=fnder(ppr1,1);
ppr1dd=fnder(ppr1,2);

[x3,y3]=buffer1d(x3,y3,n/8,1);
pp3 = csapi(x3,y3);
pp3d=fnder(pp3,1);
pp3dd=fnder(pp3,2);

xlast=.2;
yb = fnval(ppb1,nx);
yb(nx>bx1(end))=by1(end);
yt = fnval(pp3,nx);

ny=(yt+yb)/2;
mask=abs(yb-yt)>1e-7;
indexes=find(mask);
for i=indexes
    sx=nx(i);
    fun=@(yv)func(yv,sx,theta,ppb1,ppb1d,ppb1dd,bx1,by1,...
        ppr1,ppr1d,ppr1dd,rx1,ry1,pp3,pp3d,pp3dd,x3,y3,xlast,xrotate,yrotate,ST);
    
    %yv = bisectionserach(fun, yb(i), yt(i));
    yv = newtonsearch(fun, yb(i), yt(i));
    %fval=fun(yv);
    
    
    %1-sqrt(((fvaly-fval)*invh)^2+((fvalx-fval)*invh)^2)
    %fval=fun(ny(i));
    %if fval>1e-3
    %    FLAG=1;
    %    break
    %end
    
    ny(i)=yv;
end
%plot(nx,ny);
%plot(nx,yb);hold on
%plot(nx,yt)
% plot(x1,y1)
% YS=linspace(yb(i)-.2,yt(i),100);
% val=zeros(1,100);
% for j=1:100
%     val(j)=fun(YS(j));
% end
%plot(YS,val);
end

function [x0]=argmin(sx,sy,x,y)
    dist=(x-sx).^2+(y-sy).^2;
    [~,I]=min(dist);
    x0=x(I);
end

function [val]=func(yv,xv,theta,ppb,ppbd,ppbdd,bx1,by1,...
    ppr,pprd,pprdd,rx1,ry1,ppt,pptd,pptdd,x3,y3,xlast,xrotate,yrotate,ST)
    [rx,ry]=rotatepts(xv,yv,-theta,xrotate,yrotate);
    if xv<=xlast
        [xb0]=argmin(xv,yv,bx1,by1);
        [dminb,x0,~]=distancetospline(xv,yv,xb0,ppb,ppbd,ppbdd);
        %[xb0]=argmin(xv+h,yv,bx1,by1);
        %[dminbx,x0]=distancetospline(xv+h,yv,xb0,ppb,ppbd,ppbdd);
        %[xb0]=argmin(xv,yv+h,bx1,by1);
        %[dminby,x0]=distancetospline(xv,yv+h,xb0,ppb,ppbd,ppbdd);
        %ddist=1-sqrt(((dminbx-dminb)*invh)^2+((dminby-dminb)*invh)^2);
        if x0>bx1(end) || x0<bx1(1)
            'warning bot'
            n=numel(bx1)
            bx1(1)
            x0
            bx1(end)
        end
    else
        dminb=inf;
    end
    if xv>=xlast
        [xr0]=argmin(rx,ry,rx1,ry1);
        [dminr,x0,~]=distancetospline(rx,ry,xr0,ppr,pprd,pprdd);
        %[xr0]=argmin(rx+h,ry,rx1,ry1);
        %[dminrx,x0]=distancetospline(rx+h,ry,xr0,ppr,pprd,pprdd);
        %[xr0]=argmin(rx,ry+h,rx1,ry1);
        %[dminry,x0]=distancetospline(rx,ry+h,xr0,ppr,pprd,pprdd);
        %ddist=1-sqrt(((dminrx-dminr)*invh)^2+((dminry-dminr)*invh)^2);
        if x0<rx1(end) || x0>rx1(1)
            n=numel(bx1)
            'warning rot'
            rx1(end)
            x0
            rx1(1)
        end
    else
        dminr=inf;
    end
    [xt0]=argmin(xv,yv,x3,y3);
    [dmint,x0,~]=distancetospline(xv,yv,xt0,ppt,pptd,pptdd);
    %[xt0]=argmin(xv+h,yv,x3,y3);
    %[dmintx,x0]=distancetospline(xv+h,yv,xt0,ppt,pptd,pptdd);
    %[xt0]=argmin(xv,yv+h,x3,y3);
    %[dminty,x0]=distancetospline(xv,yv+h,xt0,ppt,pptd,pptdd);
    %ddist=max(abs(1-sqrt(((dmintx-dmint)*invh)^2+((dminty-dmint)*invh)^2)),abs(ddist));
    if x0>x3(end) || x0<x3(1)
        n=numel(bx1)
        'warning top'
        x3(1)
        x0
        x3(end)
    end
    phi1=min(dminb,dminr);
    phi3=dmint;
    val=(ST(1)*phi1-ST(3)*phi3);%+max(yv-ytop,0)+max(ybot-yv,0);
end


function [newx,newy]=buffer1d(x,y,pad,both)
n=numel(x);
if both
    newx=zeros(1,n+2*pad);
    newy=zeros(1,n+2*pad);
else
    newx=zeros(1,n+pad);
    newy=zeros(1,n+pad);
end
newx(pad+1:n+pad)=x;
newy(pad+1:n+pad)=y;
for i=1:pad
    newx(pad+1-i)=2*newx(pad+1)-newx(pad+1+i);
    newy(pad+1-i)=newy(pad+1+i);
    if both
        newx(n+pad+i)=2*newx(n+pad)-newx(n+pad-i);
        newy(n+pad+i)=newy(n+pad-i);
    end
end

end

function p = newtonsearch(fun, a, b)%is really the secant method.

    tol=1e-8;
    maxiter=1000;
    h=1e-8;
    invh=1e8;
    x0 = (a + b)/2;
    me=1e-11;
    x1=x0+h;
    f0=fun(x0);
    %disterror=0;
    for i=1:maxiter
        f1=fun(x1);
        %disterror=max(disterror,ddist);
        if abs(f1)<tol
            break
        end
        
        %disterror=max(disterror,ddist);
        df=(f1-f0)/(x1-x0+me);
        x0=x1;
        f0=f1;
        x1=x1-f1/(df+me);
    end
    p= x1;
end

function p = newtonsearchold(fun, a, b)%one in the paper

    tol=1e-8;
    maxiter=1000;
    h=1e-8;
    invh=1e8;
    x = (a + b)/2;
    me=1e-11;
    %disterror=0;
    for i=1:maxiter
        f=fun(x);
        %disterror=max(disterror,ddist);
        if abs(f)<tol
            break
        end
        fh=fun(x+h);
        %disterror=max(disterror,ddist);
        df=(fh-f)*invh;
        x=x-f/(df+me);

    end
    p= x;
end