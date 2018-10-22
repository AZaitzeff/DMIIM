function [lsv]=dictmapreconstcenter(point,points,phimap,num,ratio,dist)
point=point/ratio;
[z,~]=size(phimap);
[dmin,ind]=min(sum((phimap'-point).^2));
rem=mod(ind,num);


if rem==0 || rem==1 || ind<num || ind>z-num

    lsv=Inf;
    
else
    x0=points(ind,1);
    y0=points(ind,2);
    [x,y]=dist(x0,y0,point);
    lsv=sqrt(x^2+y^2);
end
end
