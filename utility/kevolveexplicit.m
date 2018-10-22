function [x,y] = kevolveexplicit(dt,nt,x,y,S,M,closed,dir,rpar)

% Takes nt steps of curvature motion.
% dt<dx^2/16 numerically

% S & M are the surface tension and mobility.
dteff=dt*S*M;
tol=1e-7;%dteff;

x=padarray(x,[0,1]);
y=padarray(y,[0,1]);
n = size(x,2);
[x,y]=bndcon(x,y,closed,n);

% Compute the normal:
for t=1:nt
    

    h=sqrt(((x(3:n) - x(1:n-2))).^2+((y(3:n) - y(1:n-2))).^2);
    k=1./sqrt((x(2:n) - x(1:n-1)).^2+(y(2:n) - y(1:n-1)).^2);
    z=2./h;
    
    x(2:n-1)=x(2:n-1)+dteff*(-x(2:n-1).*z.*(k(1:end-1)+k(2:end))+x(1:n-2).*z(1:end).*(k(1:end-1))+...
        x(3:n).*z(1:end).*(k(2:end)));
    y(2:n-1)=y(2:n-1)+dteff*(-y(2:n-1).*z.*(k(1:end-1)+k(2:end))+y(1:n-2).*z(1:end).*(k(1:end-1))+...
        y(3:n).*z(1:end).*(k(2:end)));

    [x,y]=bndcon(x,y,closed,n);


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