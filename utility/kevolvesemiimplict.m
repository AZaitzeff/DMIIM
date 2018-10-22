function [x,y] = kevolvesemiimplict(dt,nt,x,y,S,M,closed,dir,rpar)

% Takes nt steps of curvature motion.
% dt<dx^2/16 numerically

% S & M are the surface tension and mobility.
dteff=dt*S*M;

x=padarray(x,[0,1]);
y=padarray(y,[0,1]);
n = size(x,2);
[x,y]=bndcon(x,y,closed,n);



% Compute the normal:
for t=1:nt
    
    k=1./sqrt((x(2:n) - x(1:n-1)).^2+(y(2:n) - y(1:n-1)).^2);
    z=2./sqrt(((x(3:n) - x(1:n-2))).^2+((y(3:n) - y(1:n-2))).^2);
    B=zeros(n-2,3);
    B(:,2)=1+dteff*z.*(k(1:end-1)+k(2:end));
    B(1:n-3,1)=-dteff*z(2:end).*(k(2:end-1));
    B(2:n-2,3)=-dteff*z(1:end-1).*(k(2:end-1));
    Bx=B;
    By=B;
    if closed==0%neumann boundary condition derivative zero
        Bx(1,2)=Bx(1,2)-dteff*z(1)*2*(k(1));
        Bx(2,3)=Bx(2,3)+dteff*z(1)*(k(1));
        By(2,3)=By(2,3)-dteff*z(1)*k(1);

        Bx(n-2,2)=Bx(n-2,2)-dteff*z(end)*2*(k(end));
        Bx(n-3,1)=Bx(n-3,1)+dteff*z(end)*(k(end));
        By(n-3,1)=By(n-3,1)-dteff*z(end)*k(end);
    elseif closed==4
        Bx(1,2)=Bx(1,2)-dteff*z(1)*2*(k(1));
        Bx(2,3)=Bx(2,3)+dteff*z(1)*(k(1));
        By(2,3)=By(2,3)-dteff*z(1)*k(1);

        By(n-2,2)=By(n-2,2)-dteff*z(end)*2*(k(end));
        By(n-3,1)=By(n-3,1)+dteff*z(end)*(k(end));
        Bx(n-3,1)=Bx(n-3,1)-dteff*z(end)*k(end);
    end

    if dir==-1
        By(n-2,2)=1;
        By(n-3,1)=0;
    end
    if dir==1
        By(1,2)=1;
        By(2,3)=0;
    end
    
    Ax = spdiags(Bx,[-1,0,1],n-2,n-2);
    Ay = spdiags(By,[-1,0,1],n-2,n-2);
    %full(Ax)
    %full(Ay)
    x(2:n-1) = mldivide(Ax,x(2:n-1)');
    y(2:n-1) = mldivide(Ay,y(2:n-1)');
    if any(isnan([x,y])) || any(abs([x,y])>1e2)
        break
    end
        
   normpts=min(diff(x(2:n-1)).^2+ diff(y(2:n-1)).^2);
 
   if normpts<rpar 
 
      [xn,yn] = reparam(x(2:n-1),y(2:n-1));
 
      x(2:n-1) = xn; y(2:n-1) = yn;
 
   end
   [x,y]=bndcon(x,y,closed,n);

end



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