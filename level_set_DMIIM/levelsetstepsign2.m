function [phi]=levelsetstepsign2(phi,gam,h,dt,ren,eps)
    [z]=size(phi,2);
    [n,m]=size(phi{1});
    me=.000001;
    for i =1:z
        u=reshape(phi{i},[n,m])-eps;
        for t=1:ren
            ex=padarray(u,[1,0],'circular','post');
            fx=diff(ex,1,1)/h;
            ey=padarray(u,[0,1],'circular','post');
            fy=diff(ey,1,2)/h;
            
            ex=padarray(u,[1,0],'circular','pre');
            bx=diff(ex,1,1)/h;
            ey=padarray(u,[0,1],'circular','pre');
            by=diff(ey,1,2)/h;
            
            normgradc=sqrt(((fx+bx)/2).^2+((fy+by)/2).^2+me);

            
            normgrad=sqrt((fx).^2+(fy).^2+me);
            ex=padarray(fx./normgrad,[1,0],'circular','pre');
            bx=diff(ex,1,1)/h;
            ey=padarray(fy./normgrad,[0,1],'circular','pre');
            by=diff(ey,1,2)/h;

            curvature=bx+by;
            u=u+dt*gam(i)*curvature.*normgradc;
            
            %if (mod(t,100)==0)&&t~=ren
            %   temp=phi{i};
            %   temp(abs(temp)<1e-6)=0;
            %   [temp,~,~] = redistz(temp,width,flag,1/n,1/n);
            %   phi{i}=reshape(temp,[n n]);
            %end
        end
        phi{i}=u;
    end
        
    
end