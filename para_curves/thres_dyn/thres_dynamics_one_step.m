function [nx,ny,FLAG]=thres_dynamics_one_step(nx,ny,DT,ST,boost)
FLAG=0;
n=size(nx,2);
x1=[-flip(nx) nx(2:end)];
y1=[flip(ny) ny(2:end)];
x2=[-flip(nx) nx(2:end)]+.5;
y2=[flip(ny) ny(2:end)];

x3=[-flip(nx) nx(2:end-1) -flip(nx)+.5];
y3=[flip(ny) ny(2:end-1) flip(ny)];
pp1 = csapi(x1,y1);
pp2 = csapi(x2,y2);
pp3 = csapi(x3,y3);
parfor i=1:n
    sx=nx(i);
    if boost
        xs1=linspace(max((sx-8*sqrt(DT)),x1(1)),min((8*sqrt(DT)+sx),x1(end)),boost);
        ys1 = fnval(pp1,xs1);
        xs2=linspace(x2(1),min((8*sqrt(DT)+.25),x2(end)),boost);
        ys2 = fnval(pp2,xs2);
        xs3=linspace(max((sx-8*sqrt(DT)),x3(1)),min((8*sqrt(DT)+sx),x3(end)),boost);
        ys3 = fnval(pp3,xs3);
        fun=@(yv)func(yv,sx,xs1,ys1,xs2,ys2,xs3,ys3,DT,ST);
    else
        v1=(sx-8*sqrt(DT))<x1 & (8*sqrt(DT)+sx)>x1;
        v2=(8*sqrt(DT)+.25)>x2;
        v3=(sx-8*sqrt(DT))<x3 & (8*sqrt(DT)+sx)>x3;
        fun=@(yv)func(yv,sx,x1(v1),y1(v1),x2(v2),y2(v2),x3(v3),y3(v3),DT,ST);
    end
    
    yv = newtonsearch(fun, ny(i));
    %fval=fun(yv);
    
    %if fval>1e-3
    %    FLAG=1;
    %    break
    %end
    
    ny(i)=yv;
end

end


function [val]=func(yv,xv,x1,y1,x2,y2,x3,y3,DT,ST)
    
    dist1=convguass(x1,y1,xv,yv,DT,1);
    dist2=convguass(x2,y2,xv,yv,DT,1);
    dist3=convguass(x3,y3,xv,yv,DT,-1);
    
    phi1=ST(1,2)*dist2+ST(1,3)*dist3;
    phi3=ST(3,2)*dist2+ST(3,1)*dist1;
    val=(phi1-phi3);
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

function p = newtonsearch(fun, x0)

    tol=1e-8;
    maxiter=1000;
    h=1e-8;
    invh=1e8;
    x = x0;
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