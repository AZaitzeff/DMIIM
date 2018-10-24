function [lsv,phase]=dictmapreconst(point,access,pointscell,phimapcell,num,slopes,ratio,dist)
point=point/ratio;
points=pointscell{access};
phimap=phimapcell{access};
[z,~]=size(phimap);
[dmin,ind]=min(sum((phimap'-point).^2));
rem=mod(ind,num);
ind1=ind;
ind2=ind;
if rem==0 || rem==1 || ind<num || ind>z-num
    [lsv,phase]=min(point);
    point(phase)=[];
    lsv=abs(lsv-min(point));
else
    
    
    x0=points(ind,1);
    y0=points(ind,2);
    [x,y]=dist(x0,y0,point);
    slope2=slopes(1);
    slope3=slopes(2);
    ratio2=slopes(3);
    ratio3=slopes(4);
    if access==1
        if y>slope2*x
            phase=1;
        else
            phase=2;
        end
        dmin=abs(y-x*slope2)*ratio2;
    elseif access==2
        if y>slope3*x
            phase=1;
        else
            phase=3;
        end
        dmin=abs(y-x*slope3)*ratio3;
    else
        if x<=0
            phase=2;
        else
            phase=3;
        end
        dmin=abs(x);
    end
    lsv=dmin; 
end
end
