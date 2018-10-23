function [nx,ny,FLAG]=epsilon_x_spline(x,y,nx,eps,up)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds Epsilon level set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   x,y = Interface of region 1
%   nx =  x points where you want to find epsilon level set
%   eps = epsilon for epsilon level set
%   up = epsilon level set above (up=1) or below (up=-1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
%   nx,ny = Epsilon level set
%   Flag = if Flag is 1 an error has occurred

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


FLAG=0;
n=numel(nx);


[bx,by]=buffer1d(x,y,floor(n/8),1);
ppb = csapi(bx,by);
ppbd=fnder(ppb,1);
ppbdd=fnder(ppb,2);
ny = fnval(ppb,nx)+eps*up;
for i=1:n
    sx=nx(i);
    fun=@(yv)func(yv,sx,ppb,ppbd,ppbdd,bx,by,eps);
    
    yv = newtonsearch(fun, ny(i));
    fval=fun(yv);
    
    if fval>1e-3
        FLAG=1;
        break
    end
    
    ny(i)=yv;
end
end

function [x0]=argmin(sx,sy,x,y)
    dist=(x-sx).^2+(y-sy).^2;
    [~,I]=min(dist);
    x0=x(I);
end

function [val]=func(yv,xv,ppb,ppbd,ppbdd,bx,by,eps)
    [xb0]=argmin(xv,yv,bx,by);
    [dmin,x0]=distancetospline(xv,yv,xb0,ppb,ppbd,ppbdd);

    if x0>bx(end) || x0<bx(1)
        'warning bottom'
        n=numel(bx)
        bx(1)
        x0
        bx(end)
    end
    val=(dmin-eps);
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
        newy(n+pad+i)=-(newy(n+pad)-newy(n-1+pad))/(newx(n+pad)-newx(n-1+pad))*(newx(n+pad-i)-newx(n+pad))+newy(n+pad);
    end
end

end

function p = newtonsearch(fun, x0)%is really the secant method.

    tol=1e-8;
    maxiter=1000;
    h=1e-8;
    invh=1e8;
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

function p = newtonsearchold(fun,x)%method in the paper

    tol=1e-8;
    maxiter=1000;
    h=1e-8;
    invh=1e8;
    me=1e-11;
    for i=1:maxiter
        f=fun(x);
        if abs(f)<tol
            break
        end
        fh=fun(x+h);
        df=(fh-f)*invh;
        x=x-f/(df+me);

    end
    p= x;
end