function val=trap_per(Y,x)
    N=numel(Y);
    val=0;
    for i=1:N-1

        t=i+1;

        if sign(Y(i))==sign(Y(t))
            val=val+abs((Y(i)+Y(t))/2)*abs(x(t)-x(i));
        else
            total=abs(Y(i))+abs(Y(t));
            val=val+(abs(Y(i))*abs(Y(i))+abs(Y(t))*abs(Y(t)))/total*abs(x(t)-x(i))/2;
        end  
    end
end