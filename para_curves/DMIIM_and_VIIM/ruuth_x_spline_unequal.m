function [nx1,ny1,nx2,ny2,nx3,ny3]=ruuth_x_spline_unequal(x1,y1,x2,y2,x3,y3,nx1,nx2,xrotate,yrotate,angle1,angle2,ruuthfunc,ruuthfunccenter,ratio)
h=1e-8;
invh=1e8;
tol=1e-8;
tolb=1e-6;
FLAG=0;
N1=numel(nx1);
N2=numel(nx2);
angle3=(360-angle1-angle2);
%[points,phimap,xvalid]=makedictionary(angle,100,DT,ST,1e-4);
[bx1,by1]=buffer1d(x1,y1,8,1);
pp1 = csapi(bx1,by1);
pp1d=fnder(pp1,1);
pp1dd=fnder(pp1,2);
%theta=2*pi/3;

theta2=(180-angle2/2)/180*pi;
theta3=(180-angle3/2)/180*pi;


[rx2,ry2]=rotatepts(x2,y2,-theta2,xrotate,yrotate);
ppr2 = csapi(rx2,ry2);
ppr2d=fnder(ppr2,1);
ppr2dd=fnder(ppr2,2);

[bx2,by2]=buffer1d(x2(1:N1),y2(1:N1),8,0);
pp2 = csapi(bx2,by2);
pp2d=fnder(pp2,1);
pp2dd=fnder(pp2,2);


[rx3,ry3]=rotatepts(x3,y3,theta3,xrotate,yrotate);
ppr3 = csapi(rx3,ry3);
ppr3d=fnder(ppr3,1);
ppr3dd=fnder(ppr3,2);

[bx3,by3]=buffer1d(x3(1:N2),y3(1:N2),8,0);
pp3 = csapi(bx3,by3);
pp3d=fnder(pp3,1);
pp3dd=fnder(pp3,2);
xlast2=xrotate/2;
xlast3=.5-(.5-xrotate)/2;
dist=@(x,y) phimapfunc(x,y,theta2,theta3,pp1,pp1d,pp1dd,x1,y1,...
    ppr2,ppr2d,ppr2dd,rx2,ry2,pp2,pp2d,pp2dd,bx2,by2,xlast2,...
    ppr3,ppr3d,ppr3dd,rx3,ry3,pp3,pp3d,pp3dd,bx3,by3,xlast3,xrotate,yrotate);


ux2=flip(x2(N1:end));
uy2=flip(y2(N1:end));

[uy2,ux2]=buffer1d(uy2,ux2,8,0);

ppu2 = csapi(uy2,ux2);
ppu2d=fnder(ppu2,1);
ppu2dd=fnder(ppu2,2);

ux3=flip(x3(N2:end));
uy3=flip(y3(N2:end));

[uy3,ux3]=buffer1d(uy3,ux3,8,0);

ppu3 = csapi(uy3,ux3);
ppu3d=fnder(ppu3,1);
ppu3dd=fnder(ppu3,2);

ylast2=-.5-(-.5-y2(N1))/2;
ylast3=-.5-(-.5-y3(N2))/2;

distside=@(x,y) phimapfuncside(x,y,theta2,theta3,pp1,pp1d,pp1dd,x1,y1,...
    ppr2,ppr2d,ppr2dd,rx2,ry2,ppu2,ppu2d,ppu2dd,ux2,uy2,ylast2,...
    ppr3,ppr3d,ppr3dd,rx3,ry3,ppu3,ppu3d,ppu3dd,ux3,uy3,ylast3,xrotate,yrotate);

center=@(x,y)ruuthfunccenter(dist(x,y));

[xc,yc] = gridsearchsearch2d(center, xrotate, yrotate);
nx1=nx1*xc/nx1(end);
nx2=.5-(.5-nx2)*(.5-xc)/(.5-nx2(end));

yt = fnval(pp1,nx1);
yb = fnval(pp2,nx1);
yb(nx1>bx2(end))=by2(end);

ny1=(yt+yb)/2;
mask=abs(yb-yt)>tolb;
%indexes=find(mask);
dictmap=@(x,y)ruuthfunc(dist(x,y),1);
Nz=numel(ny1);
for i=1:Nz-1
    if mask(i)
    sx=nx1(i);
    f=@(y) redistfunc(sx,y,dictmap,1);
    yv=newtonsearch(f,ny1(i),h,invh,tol,yb(i),yt(i)); 
    ny1(i)=yv;
    end
end
ny1(end)=yc;
yt = fnval(pp1,nx2);
yb = fnval(pp3,nx2);
yb(nx1>bx3(end))=by3(end);

ny2=(yt+yb)/2;
mask=abs(yb-yt)>tolb;
%indexes=find(mask);
dictmap=@(x,y)ruuthfunc(dist(x,y),2);
Nz=numel(ny2);
for i=1:Nz-1
    if mask(i)
    sx=nx2(i);
    f=@(y) redistfunc(sx,y,dictmap,1);
    %ys=linspace(ny2(i)-.05,ny2(i)+.05,501);
    %master=zeros(1,501);

    %for j=1:501    
    %    master(j)=f(ys(j));
    %end
    %plot(ys,master);hold on;plot(ys(251),master(251),'*')
    %yv=newtonsearch(f,ny2(i),h,invh,tol,yb(i),yt(i));
    %return
    %plot(x1,y1);hold on;plot(sx,ny1(i),'*');
    
   
    yv=newtonsearch(f,ny2(i),h,invh,tol,yb(i),yt(i));
    ny2(i)=yv;
    end
end
ny2(end)=yc;
avgnorm=(sqrt((ny2(end)-ny2(end-1))^2+(nx2(end)-nx2(end-1))^2)+...
    sqrt((ny1(end)-ny1(end-1))^2+(nx1(end)-nx1(end-1))^2))/2;
%avgnorm=2.7151e-04;
border=y2(end);
M=ceil((yc-border)/avgnorm);
ny3=linspace(yc,border,M);
xb = fnval(ppu2,ny3);
xb(ny3>uy2(end))=ux2(end);
xt = fnval(ppu3,ny3);
xt(ny3>uy3(end))=ux3(end);


nx3=(xt+xb)/2;
dictmap=@(x,y)ruuthfunc(distside(x,y),3);
mask=abs(xb-xt)>tolb;
%indexes=find(mask);
nx3(1)=xc;
Nz=numel(nx3);
for i=2:Nz
    if mask(i)
    sy=ny3(i);
    f=@(x) redistfunc(x,sy,dictmap,3);
    %xs=linspace(.25-.01,.25+.01,101);
    %master=zeros(1,101);
    %for j=1:101    
    %    master(j)=f(xs(j));
    %end
    %plot(xs,master);%hold on;plot(xs(51),master(51),'*')
    %return
    %plot(x1,y1);hold on;plot(sx,ny1(i),'*');
    xv=newtonsearch(f,nx3(i),h,invh,tol,xb(i),xt(i));
    nx3(i)=xv;
    end
end

end

function [x0]=argmin(sx,sy,x,y)
    dist=(x-sx).^2+(y-sy).^2;
    [~,I]=min(dist);
    x0=x(I);
end

function [val]=redistfunc(x,y,func,actphase)
    [lsv,phase]=func(x,y);
     val=lsv*(-1)^(phase==actphase);
end

function [phi]=phimapfunc(xc,yc,theta2,theta3,pp1,pp1d,pp1dd,x1,y1,...
    ppr2,ppr2d,ppr2dd,rx2,ry2,pp2,pp2d,pp2dd,x2,y2,xlast2,...
    ppr3,ppr3d,ppr3dd,rx3,ry3,pp3,pp3d,pp3dd,x3,y3,xlast3,xrotate,yrotate)
    phi=zeros(3,1);
    
    
    [x0]=argmin(xc,yc,x1,y1);
    [dmin1,~,y0]=distancetospline(xc,yc,x0,pp1,pp1d,pp1dd);
    if yc<y0
        dmin1=-dmin1;
    end
    
    if xc<=xlast2
        [x0]=argmin(xc,yc,x2,y2);
        [dmin2,~,y0]=distancetospline(xc,yc,x0,pp2,pp2d,pp2dd);
        if yc>y0
            dmin2=-dmin2;
        end
    else
        [xcr2,ycr2]=rotatepts(xc,yc,-theta2,xrotate,yrotate);
        [x0]=argmin(xcr2,ycr2,rx2,ry2);
        [dmin2,~,y0]=distancetospline(xcr2,ycr2,x0,ppr2,ppr2d,ppr2dd);
        if ycr2<y0
            dmin2=-dmin2;
        end
    end
    
    if xc>=xlast3
        [x0]=argmin(xc,yc,x3,y3);
        [dmin3,~,y0]=distancetospline(xc,yc,x0,pp3,pp3d,pp3dd);
        if yc>y0
            dmin3=-dmin3;
        end
    else
        [xcr3,ycr3]=rotatepts(xc,yc,theta3,xrotate,yrotate);
        [x0]=argmin(xcr3,ycr3,rx3,ry3);
        [dmin3,~,y0]=distancetospline(xcr3,ycr3,x0,ppr3,ppr3d,ppr3dd);
        if ycr3<y0
            dmin3=-dmin3;
        end
    end

    
    phi(1)=dmin1;
    phi(2)=dmin2;
    phi(3)=dmin3;
end


function [phi]=phimapfuncside(xc,yc,theta2,theta3,pp1,pp1d,pp1dd,x1,y1,...
    ppr2,ppr2d,ppr2dd,rx2,ry2,ppu2,ppu2d,ppu2dd,x2,y2,ylast2,...
    ppr3,ppr3d,ppr3dd,rx3,ry3,ppu3,ppu3d,ppu3dd,x3,y3,ylast3,xrotate,yrotate)
    phi=zeros(3,1);
    
    
    [x0]=argmin(xc,yc,x1,y1);
    [dmin1,~,y0]=distancetospline(xc,yc,x0,pp1,pp1d,pp1dd);
    if yc>y0
        dmin1=-dmin1;
    end
    
    if yc<=ylast2
        [y0]=argmin(yc,xc,y2,x2);
        [dmin2,~,x0]=distancetospline(yc,xc,y0,ppu2,ppu2d,ppu2dd);
        if xc<x0
            dmin2=-dmin2;
        end
    else
        [xcr2,ycr2]=rotatepts(xc,yc,-theta2,xrotate,yrotate);
        [x0]=argmin(xcr2,ycr2,rx2,ry2);
        [dmin2,~,y0]=distancetospline(xcr2,ycr2,x0,ppr2,ppr2d,ppr2dd);
        if ycr2>y0
            dmin2=-dmin2;
        end
    end
    
    if yc<=ylast3
        [y0]=argmin(yc,xc,y3,x3);
        [dmin3,~,x0]=distancetospline(yc,xc,y0,ppu3,ppu3d,ppu3dd);
        if xc>x0
            dmin3=-dmin3;
        end
    else
        [xcr3,ycr3]=rotatepts(xc,yc,theta3,xrotate,yrotate);
        [x0]=argmin(xcr3,ycr3,rx3,ry3);
        [dmin3,~,y0]=distancetospline(xcr3,ycr3,x0,ppr3,ppr3d,ppr3dd);
        if ycr3>y0
            dmin3=-dmin3;
        end
    end

    
    phi(1)=dmin1;
    phi(2)=dmin2;
    phi(3)=dmin3;
end

function val=ruuth_projection(y,f,phimap,points)
    phi=f(y);
    [~,I]=min(sum((phimap-phi).^2,2));
    %I = knnsearch(phimap,phi)
    %[~,I]=min(sum((phimap-phi).^2,2));
    region=points(I,3);
    if region==1
        val=-1;
    elseif region==3
        val=1;
    elseif region==0
        val=0;
    else
        error('oh no');
    end
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

function [x,y] = gridsearchsearch2d(func, x0, y0)
    x=x0;
    y=y0;
    N=21;
    fac=1;
    for z=1:6
    thex=linspace(x-.05*fac,x+.05*fac,N);
    they=linspace(y-.05*fac,y+.05*fac,N);
    master=zeros(11,11);
    for i=1:N
        for j=1:N
            master(i,j)=func(thex(i),they(j));
        end
    end
    [~,I]=min(master(:));
    [row,col]=ind2sub([N,N],I);
    x=thex(row);
    y=they(col);
    fac=fac/10;
    end
    %[X,Y]=meshgrid(thex,they);
    %surf(master);
end

function [x,y] = newtonsearch2d(fun, x0, y0)

    tol=1e-5;
    maxiter=2;
    h=1e-4;
    invh=1e4;
    x = x0;
    y = y0;
    me=1e-11;
    %disterror=0;
    for i=1:maxiter
        f=fun(x,y);
        
        %disterror=max(disterror,ddist);
        if abs(f)<tol
            break
        end
        fpx=fun(x+h,y);
        fnx=fun(x-h,y);
        fny=fun(x,y-h);
        fpy=fun(x,y+h);
        fpxpy=fun(x+h,y+h);
        fnxpy=fun(x-h,y+h);
        fnxny=fun(x-h,y-h);
        fpxny=fun(x+h,y-h);
        
        %disterror=max(disterror,ddist);
        dfx=(fpx-fnx)*invh/2;
        dfy=(fpy-fny)*invh/2;
        dfxx=(fpx-2*f+fnx)*invh^2;
        dfyy=(fpy-2*f+fny)*invh^2;
        dfxy=(fpxpy-fnxpy-fpxny+fnxny)*invh^2/4;
        dethes=dfxx*dfyy-dfxy^2+me;
        x=x-(dfx*dfyy-dfy*dfxy)/dethes;
        y=y-(-dfx*dfxy+dfy*dfxx)/dethes;

    end
end

function p = newtonsearch(fun, x0,h,invh,tol,limb,limt)%secant method

    maxiter=50;
    %h=1e-3;
    %invh=1e3;
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
        df=(f1-f0)/(x1-x0);
        f0=f1;
        x0=x1;
        x1=x1-f1/(df+me);

    end
    if i==maxiter || limb-x>1e-5 || x-limt>1e-5
        limb-x
        x-limt
        p = bisectionserach(fun, limb,limt,tol);
        warning('did bisection search')

        
    else
        p= x1;
    end
    
end

function p = newtonsearchold(fun, x0,h,invh,tol,limb,limt)%as in paper

    maxiter=50;
    %h=1e-3;
    %invh=1e3;
    x=x0;
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
    if i==maxiter || limb-x>1e-5 || x-limt>1e-5
        limb-x
        x-limt
        p = bisectionserach(fun, limb,limt,tol);
        warning('did bisection search')

        
    else
        p= x;
    end
    
end

function p = bisectionserach(f, a, b,tol)
    maxiter=2000;
    for i=1:maxiter
        c = (a + b)/2;
        val=f(c);
        if abs(val)<tol || (b-a)<1e-9
            abs(val)
            break
        end
        if val>0
           a=c;
        else
           b=c; 
        end

    end
    p= c;
end