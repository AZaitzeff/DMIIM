function [phi]=levelsetstepsign3(phi,gam,h,dt,ren,width,flag)
    [~,z]=size(phi);
    [n,~,~]=size(phi{1});
    me=1e-7;
    for i =1:z

        for t=1:ren
            ex=padarray(phi{i},[1,0,0],'circular','post');
            fx=diff(ex,1,1)/h;
            ey=padarray(phi{i},[0,1,0],'circular','post');
            fy=diff(ey,1,2)/h;
            ez=padarray(phi{i},[0,0,1],'circular','post');
            fz=diff(ez,1,3)/h;
            
            ex=padarray(phi{i},[1,0,0],'circular','pre');
            bx=diff(ex,1,1)/h;
            ey=padarray(phi{i},[0,1,0],'circular','pre');
            by=diff(ey,1,2)/h;
            ez=padarray(phi{i},[0,0,1],'circular','pre');
            bz=diff(ez,1,3)/h;
            
            normgradc=sqrt(((fx+bx)/2).^2+((fy+by)/2).^2+((fz+bz)/2).^2+me);

            
            normgrad=sqrt((fx).^2+(fy).^2+(fz).^2+me);
            ex=padarray(fx./normgrad,[1,0,0],'circular','pre');
            bx=diff(ex,1,1)/h;
            ey=padarray(fy./normgrad,[0,1,0],'circular','pre');
            by=diff(ey,1,2)/h;
            ez=padarray(fz./normgrad,[0,0,1],'circular','pre');
            bz=diff(ez,1,3)/h;
            
            curvature=bx+by+bz;
            phi{i}=phi{i}+dt*gam(i)*curvature.*normgradc;
            if (mod(t,400)==0)&&t~=ren
               temp=phi{i};
               temp(abs(temp)<1e-6)=0;
               [temp,~,~] = redistz3(temp,width,flag,n,1/n,1/n,1/n);
               phi{i}=reshape(temp,[n n n]);
            end
       
        end
    end
        
    
end