function [x,y] = kevolvecn(dt,nt,x,y,S,M,closed,dir,rpar)

% Takes nt steps of curvature motion.
% dt<dx^2/16 numerically

% S & M are the surface tension and mobility.
dteff=dt*S*M;
tol=dteff;

x=padarray(x,[0,1]);
y=padarray(y,[0,1]);
n = size(x,2);
[x,y]=bndcon(x,y,closed,n);
MAXITER=400;
%iters=[];
% Compute the normal:
for t=1:nt
    xl=x;
    yl=y;
    Lp=0;
    h=sqrt(((x(3:n) - x(1:n-2))).^2+((y(3:n) - y(1:n-2))).^2);
    k=1./sqrt((x(2:n) - x(1:n-1)).^2+(y(2:n) - y(1:n-1)).^2);
    z=2./h;
    for iter=1:MAXITER
        hl=sqrt(((xl(3:n) - xl(1:n-2))).^2+((yl(3:n) - yl(1:n-2))).^2);
        L=sum(hl);
        
        
        if abs(L-Lp)<tol
            break
        end
        Lp=L;
        kl=1./sqrt((xl(2:n) - xl(1:n-1)).^2+(yl(2:n) - yl(1:n-1)).^2);
        %[v1,I1]=min(hl)
        %[v2,I2]=min(sqrt((xl(2:n) - xl(1:n-1)).^2+(yl(2:n) - yl(1:n-1)).^2))
        zl=2./(hl);
        B=zeros(n-2,3);
        B(:,2)=1+dteff*zl.*(kl(1:end-1)+kl(2:end))/2;
        B(1:n-3,1)=-dteff*zl(2:end).*(kl(2:end-1))/2;
        B(2:n-2,3)=-dteff*zl(1:end-1).*(kl(2:end-1))/2;
        %Bx=B;
        %By=B;
        
        A = spdiags(B,[-1,0,1],n-2,n-2);
        Ax=A;
        Ay=A;
        if closed==0%neumann boundary condition derivative zero
            Ax(1,1)=Ax(1,1)-dteff*zl(1)*(kl(1));
            Ax(1,2)=Ax(1,2)+dteff*zl(1)*(kl(1))/2;
            Ay(1,2)=Ay(1,2)-dteff*zl(1)*kl(1)/2;

            Ax(n-2,n-2)=Ax(n-2,n-2)-dteff*zl(end)*(kl(end));
            Ax(n-2,n-3)=Ax(n-2,n-3)+dteff*zl(end)*(kl(end))/2;
            Ay(n-2,n-3)=Ay(n-2,n-3)-dteff*zl(end)*kl(end)/2;
            
            %Bx(1,2)=Bx(1,2)-dteff*zl(1)*(kl(1));
            %Bx(2,3)=Bx(2,3)+dteff*zl(1)*(kl(1))/2;
            %By(2,3)=By(2,3)-dteff*zl(1)*kl(1)/2;

            %Bx(n-2,2)=Bx(n-2,2)-dteff*zl(end)*(kl(end));
            %Bx(n-3,1)=Bx(n-3,1)+dteff*zl(end)*(kl(end))/2;
            %By(n-3,1)=By(n-3,1)-dteff*zl(end)*kl(end)/2;
        elseif closed==4
            Ax(1,1)=Ax(1,1)-dteff*zl(1)*(kl(1));
            Ax(1,2)=Ax(1,2)+dteff*zl(1)*(kl(1))/2;
            Ay(1,2)=Ay(1,2)-dteff*zl(1)*kl(1)/2;

            Ay(n-2,n-2)=Ay(n-2,n-2)-dteff*zl(end)*(kl(end));
            Ay(n-2,n-3)=Ay(n-2,n-3)+dteff*zl(end)*(kl(end))/2;
            Ax(n-2,n-3)=Ax(n-2,n-3)-dteff*zl(end)*kl(end)/2;

            %By(n-2,2)=By(n-2,2)-dteff*zl(end)*(kl(end));
            %By(n-3,1)=By(n-3,1)+dteff*zl(end)*(kl(end))/2;
            %Bx(n-3,1)=Bx(n-3,1)-dteff*zl(end)*kl(end)/2;
        end
        %full(Ay)
        xtemp=x(2:n-1)+dteff/2*(-x(2:n-1).*z.*(k(1:end-1)+k(2:end))+x(1:n-2).*z(1:end).*(k(1:end-1))+...
            x(3:n).*z(1:end).*(k(2:end)));
        ytemp=y(2:n-1)+dteff/2*(-y(2:n-1).*z.*(k(1:end-1)+k(2:end))+y(1:n-2).*z(1:end).*(k(1:end-1))+...
            y(3:n).*z(1:end).*(k(2:end)));
        xl(2:n-1) = mldivide(Ax,xtemp');
        yl(2:n-1) = mldivide(Ay,ytemp');
        [xl,yl]=bndcon(xl,yl,closed,n);
        %plot(xl,yl);hold on
    end
    %iters(t)=iter;
    if iter==MAXITER
        hl=sqrt(((xl(3:n) - xl(1:n-2))).^2+((yl(3:n) - yl(1:n-2))).^2);
        L=sum(hl);
        abs(L-Lp)
        iter
    end
    %'here'
    x=xl;
    y=yl;
    if any(isnan([x,y])) || any(abs([x,y])>1e2)
        break
    end
        
   normpts=min(diff(x(2:n-1)).^2+ diff(y(2:n-1)).^2);
   if normpts<rpar 
 
      [xn,yn] = reparam(x(2:n-1),y(2:n-1));
 
     x(2:n-1) = xn; y(2:n-1) = yn;
 
   end
   %[x,y]=bndcon(x,y,closed,n);

end


%plot(1:nt,iters);
x=x(2:n-1);
y=y(2:n-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x,y]=bndcon(x,y,closed,n)
    if closed==0%neumann boundary condition derivative zero
        x(1)=2*x(2)-x(3);
        y(1)=y(3);
        x(n)=2*x(n-1)-x(n-2);
        y(n)=y(n-2);
    elseif closed==2  %linear
        x(1)=2*x(2)-x(3);
        y(1)=2*y(2)-y(3);
        x(n)=2*x(n-1)-x(n-2);
        y(n)=2*y(n-1)-y(n-2);
    elseif closed==4
        x(1)=2*x(2)-x(3);
        y(1)=y(3);
        x(n)=x(n-2);
        y(n)=2*y(n-1)-y(n-2);
    elseif closed==6
        x(1)=2*x(2)-x(3);
        y(1)=y(3);
        x(n)=2*x(n-1)-x(n-2);
        y(n)=2*y(n-1)-y(n-2);
    else %periodic
        x(1)=x(n-1);
        y(1)=y(n-1);
        x(n)=x(2);
        y(n)=y(2);
    end
end