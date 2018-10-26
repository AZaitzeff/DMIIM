function [lsv]=dictmaprecgrid(point,points,phimap,num,numpts,slopes,ratio,dist)
point=point/ratio;
[~,ind]=min(sum((phimap'-point).^2));
%rem=mod(ind,num);
rem2=mod(ind,numpts);
if rem2>numpts-num %|| max(abs(point))>.3
    [dmax,phase]=max(point);
    lsv=point-dmax;
    point(phase)=[];
    lsv(phase)=(dmax-max(point));
else
    lsv=zeros(1,3);
    
    
    x0=points(ind,1);
    y0=points(ind,2);
    [x,y]=dist(x0,y0,point);
    slope2=slopes(1);
    slope3=slopes(2);
    ratio2=slopes(3);
    ratio3=slopes(4);
        
    if y>-x/slope2
        dmin2=abs(y-x*slope2)*ratio2;
    else
        dmin2=sqrt(x^2+y^2);
    end
    if y>-x/slope3
        dmin3=abs(y-x*slope3)*ratio3;
    else
        dmin3=sqrt(x^2+y^2);
    end
    if y<0
        dminx=abs(x);
    else
        dminx=sqrt(x^2+y^2);
    end
   lsv(1)=-min(dmin2,dmin3);
   lsv(2)=-min(dmin2,dminx);
   lsv(3)=-min(dminx,dmin3);
        
   if x<=0
        if y>slope2*x
            
            phase=1;
        else
            phase=2;
        end
        
   else
       if y>slope3*x
            phase=1;
        else
            phase=3;
        end
   end
   lsv(phase)=-lsv(phase);

end
end