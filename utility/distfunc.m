function [x,y]=distfunc(x0,y0,point,theta1,theta2,theta3,dist1,dist2,dist3)

x=x0;
y=y0;

MAXITER=59;
tol=1e-10;
temp=inf;
bestx=0;
besty=0;
me=1e-11;
%---------------
%this code tested formulas in sec 6.2 vs finite differences
%  h=1e-4;
%  invh=1e4;
%  [xr,yr]=rotatepts(x,y,-theta3,0,0);
%  [hval,hx,hy,hxx,hxy,hyy]=dist3(xr,yr,-theta3);
% 
% [xr,yr]=rotatepts(x,y+h,-theta3,0,0);
%  [hfy,~,~,~,~,~]=dist3(xr,yr,-theta3);
%  
% [xr,yr]=rotatepts(x,y-h,-theta3,0,0);
%  [hby,~,~,~,~,~]=dist3(xr,yr,-theta3);
%  
% [xr,yr]=rotatepts(x+h,y,-theta3,0,0);
%  [hfx,~,~,~,~,~]=dist3(xr,yr,-theta3);
%  
% [xr,yr]=rotatepts(x-h,y,-theta3,0,0);
%  [hbx,~,~,~,~,~]=dist3(xr,yr,-theta3);
%  
%  [xr,yr]=rotatepts(x-h,y-h,-theta3,0,0);
%  [hbxby,~,~,~,~,~]=dist3(xr,yr,-theta3);
%  
%  [xr,yr]=rotatepts(x-h,y+h,-theta3,0,0);
%  [hbxfy,~,~,~,~,~]=dist3(xr,yr,-theta3);
%  
%   [xr,yr]=rotatepts(x+h,y-h,-theta3,0,0);
%  [hfxby,~,~,~,~,~]=dist3(xr,yr,-theta3);
%  
%  [xr,yr]=rotatepts(x+h,y+h,-theta3,0,0);
%  [hfxfy,~,~,~,~,~]=dist3(xr,yr,-theta3);
% % 
% % 
% %
% 
% hxy-(hbxby-hfxby-hbxfy+hfxfy)*invh^2/4
% hxx-(hfx-2*hval+hbx)*invh^2
% hyy-(hfy-2*hval+hby)*invh^2
% 
% hx-(hfx-hbx)*invh/2
% hy-(hfy-hby)*invh/2
% return
%---------------
%vals=[];
count=1;
for i=1:MAXITER
    
        
    objx=0;
    objy=0;
    objxx=0;
    objxy=0;
    objyy=0;
    n1=0;
    nx=0;
    ny=0;
    [objx,objy,objxx,objxy,objyy,n1,nx,ny]=adder(x,y,point(1),dist1,-theta1,objx,objy,objxx,objxy,objyy,n1,nx,ny);
    [objx,objy,objxx,objxy,objyy,n1,nx,ny]=adder(x,y,point(2),dist2,-theta2,objx,objy,objxx,objxy,objyy,n1,nx,ny);
    [objx,objy,objxx,objxy,objyy,n1,nx,ny]=adder(x,y,point(3),dist3,theta3,objx,objy,objxx,objxy,objyy,n1,nx,ny);
    newval=max(abs(objx),abs(objy));
    if newval<tol
        break
    elseif newval<temp
        temp=newval;
        bestx=x;
        besty=y;
    elseif isnan(newval)
        error('newval is nan')
    end
    
    
    
    %vals(i)=max(abs(objx)/(sqrt(n1*nx)),abs(objy)/(sqrt(n1*ny)));
    %if i==1
    %    max(abs(objx)/(sqrt(n1*nx)),abs(objy)/(sqrt(n1*ny)))
    %end
    %d=[[objxx,objxy];[objxy,objyy]]\[objx;objy];
    %x=x-d(1);
    %y=y-d(2);
    if newval<1e-4 || count<10
        x=x-(objyy*objx-objxy*objy)/(objxx*objyy-objxy^2);
        y=y-(-objxy*objx+objxx*objy)/(objxx*objyy-objxy^2);
    else
        count=1;
        x=x0+randn()*.01;
        y=y0+randn()*.01;   
    end
    count=count+1;
    
end
%vals(i)=max(abs(objx)/(sqrt(n1*nx)),abs(objy)/(sqrt(n1*ny)));
if i==MAXITER
    x=bestx;
    y=besty;
    if temp>1e-3
        temp
        %save('file.mat','x0','y0','point');
    end
    
end
%max(abs(objx)/(sqrt(n1*nx)),abs(objy)/(sqrt(n1*ny)))
%semilogy(vals);
%plot(valsx(6:end),valsy(6:end));hold on
%plot(valsx(1),valsy(1),'*')
    %[xr,yr]=rotatepts(x,y,-theta1,0,0);
    %[g,gx,gy,gxx,gxy,gyy,gx0x0,gx0x,gx0y,x0x,x0y]=dist1(xr,yr);
    %g
    %[xr,yr]=rotatepts(x,y,-theta2,0,0);
    %[g,gx,gy,gxx,gxy,gyy,gx0x0,gx0x,gx0y,x0x,x0y]=dist2(xr,yr);
    %g
    %[xr,yr]=rotatepts(x,y,theta3,0,0);
    %[g,gx,gy,gxx,gxy,gyy,gx0x0,gx0x,gx0y,x0x,x0y]=dist3(xr,yr);
    %g

end

function [objx,objy,objxx,objxy,objyy,n1,nx,ny]=adder(x,y,p,dist,theta,objx,objy,objxx,objxy,objyy,n1,nx,ny)
    
    [xr,yr]=rotatepts(x,y,theta,0,0);
    [h,hx,hy,hxx,hxy,hyy]=dist(xr,yr,theta);
    
    gp=(h-p);
    objx=objx+hx*gp;
    objy=objy+hy*gp;
    objxx=objxx+hxx*gp+hx^2;
    objxy=objxy+hxy*gp+hx*hy;
    objyy=objyy+hyy*gp+hy^2;
    n1=n1+gp^2;
    nx=nx+hx^2;
    ny=ny+hy^2;
end

function [hx,hy,hxx,hxy,hyy]=distder(g,gx,gy,gx0x0,gx0x,gx0y,x0x,x0y,theta)
Rxx=cos(theta);
Rxy=-sin(theta);
Ryx=sin(theta);
Ryy=cos(theta);
hx=gx*Rxx+gy*Ryx;
hy=gx*Rxy+gy*Ryy;
hxx=(x0x*Rxx+x0y*Ryx)*(2*gx0x*Rxx+2*gx0y*Ryx+gx0x0*(x0x*Rxx+x0y*Ryx));
hyy=(x0x*Rxy+x0y*Ryy)*(2*gx0x*Rxy+2*gx0y*Ryy+gx0x0*(x0x*Rxy+x0y*Ryy));
hxy=(x0x*Rxx+x0y*Ryx)*(gx0x*Rxy+gx0y*Ryy)+...
    (x0x*Rxy+x0y*Ryy)*(gx0x*Rxx+gx0y*Ryx)+...
    gx0x0*(x0x*Rxy+x0y*Ryy)*(x0x*Rxx+x0y*Ryx);
end
