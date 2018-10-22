function [u]=upbuffer(u,lin)


[~,n]=size(u);
if lin==1
    u(2)=(2*u(3)-u(4));
    u(1)=(2*u(4)-u(5));
    u(n-1)=(2*u(n-2)-u(n-3));
    u(n)=(2*u(n-2)-u(n-4));
else
    u(1)=u(5);
    u(2)=u(4);

    u(n-1)=u(n-3);
    u(n)=u(n-4);
end