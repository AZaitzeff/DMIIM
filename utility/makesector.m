function [x,u]=makesector(angle,N)
    theta1=angle*pi/180;
    %ratio=sin(theta1/2)/sin(120*pi/180/2);
    M=cot(theta1/2);
    endpt=10;
    x=linspace(0,endpt,N);
    h=x(2)-x(1);
    u=M*exp(-(M^2+1)/4*x.^2).*(x.*exp((M^2+1)/4*x.^2).*erf(x*sqrt(M^2+1)/2)+2/(sqrt(pi)*sqrt(M^2+1)));
    tol=1e-10;
    tol2=1e-12;
    for i=1:10000
        %u0=u;
        utemp=u;
        u=buffer(u,0);
        
        du=center(u);
        du=du(3:N+2)/h;
        du(N)=M;
        f=(1+du.^2)/2;
        B=zeros(N,3);
        B(:,2)=-f-2/h^2;
        B(1:N-1,1)=1/h^2-x(2:end).*f(2:end)/(2*h);
        B(2:N,3)=1/h^2+x(1:end-1).*f(1:end-1)/(2*h);

        B(2,3)=B(2,3)+1/h^2-x(1).*f(1)/(2*h);

        B(N,2)=1;
        B(N-1,1)=0;

        A=spdiags(B,[-1,0,1],N,N);
        vec=zeros(N,1);
        vec(N)=M*x(end);

        u = mldivide(A,vec)';
        uc=buffer(u,0);
        du=center(uc)/h;
        ddu=centerxx(uc)/h^2;
        factor=max(abs(ddu(3:N+1)-(uc(3:N+1)-x(1:end-1).*du(3:N+1)).*(1+du(3:N+1).^2)/2));
        
        change=max(abs(utemp-u));
        if factor<tol || change<tol2
            break;
        end
    end
    %if factor>tol
    %    'warning'
    %end
end


