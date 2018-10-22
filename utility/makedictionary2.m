function [points,phimap,slopes,numpts,num,dist]=makedictionary2(angle1,angle2,DT,gam,stepsize,numlayer,single)
%gam = mu*sigma

angle3=(360-angle1-angle2);
slope2=-tan((angle2-90)*pi/180);
slope3=tan((angle3-90)*pi/180);

ratio2=cos((angle2-90)*pi/180);
ratio3=cos((angle3-90)*pi/180);
[odex,odey]=makesector(angle1,20000);
%[odey,~,odex] = givenm(angle1);
%odey=odey';
odex=odex(1:4:end);
odey=odey(1:4:end);
N=size(odex,2);
transx1=odex*sqrt(gam(1)*DT);
transy1=odey*sqrt(gam(1)*DT);
transx1=[-flip(transx1) transx1(2:end)];
transy1=[flip(transy1) transy1(2:end)];
theta1=(180-angle2-angle1/2)/180*pi;
%[transx3,transy3]=rotatepts(transx3,transy3,theta3,0,0);
pp1 = csapi(transx1,transy1);
pp1d=fnder(pp1,1);
pp1dd=fnder(pp1,2);


[odex,odey]=makesector(angle2,20000);
odex=odex(1:4:end);
odey=odey(1:4:end);
%[odey,~,odex] = givenm(angle2);
%odey=odey';
%odex=odex(1:100:end);
%odey=odey(1:100:end);
transx2=(odex*sqrt(gam(2)*DT));
transy2=(odey*sqrt(gam(2)*DT));
transx2=[-flip(transx2) transx2(2:end)];
transy2=[flip(transy2) transy2(2:end)];
pp2 = csapi(transx2,transy2);
pp2d=fnder(pp2,1);
pp2dd=fnder(pp2,2);

[odex,odey]=makesector(angle3,20000);
odex=odex(1:4:end);
odey=odey(1:4:end);
%[odey,~,odex] = givenm(angle3);
%odey=odey';
%odex=odex(1:100:end);
%odey=odey(1:100:end);
transx3=(odex*sqrt(gam(3)*DT));
transy3=(odey*sqrt(gam(3)*DT));
transx3=[-flip(transx3) transx3(2:end)];
transy3=[flip(transy3) transy3(2:end)];
pp3 = csapi(transx3,transy3);
pp3d=fnder(pp3,1);
pp3dd=fnder(pp3,2);



theta2=(180-angle2/2)/180*pi;

theta3=(180-angle3/2)/180*pi;

minx1=min(transx1)*.9;
minx2=min(transx2)*.9;
minx3=min(transx3)*.9;

%M=ceil(.2/stepsize);
%xvalid=.1;
%xb=linspace(-xvalid,xvalid,M);
%xvalid=xvalid*.8;
%xs,ys=meshgrid(xb,xb);
%xlast=xb(end-1);
%phimap=zeros(M,M,3);
%points=zeros(M,M,2);
slopes=zeros(1,4);
slopes(1)=slope2;
slopes(2)=slope3;
slopes(3)=ratio2;
slopes(4)=ratio3;
numpts=ceil(.22/stepsize);
slope1r=tan((90-angle1/2)*pi/180);
ratio1r=cos((90-angle1/2)*pi/180);
slope2r=tan((90-angle2/2)*pi/180);
ratio2r=cos((90-angle2/2)*pi/180);
slope3r=tan((90-angle3/2)*pi/180);
ratio3r=cos((90-angle3/2)*pi/180);

total=(numlayer*2+1)*numpts;

dist1=@(xz,yz,theta)distanceall(xz,yz,theta,transx1,transy1,pp1,pp1d,pp1dd,minx1,slope1r,ratio1r);
dist2=@(xz,yz,theta)distanceall(xz,yz,theta,transx2,transy2,pp2,pp2d,pp2dd,minx2,slope2r,ratio2r);
dist3=@(xz,yz,theta)distanceall(xz,yz,theta,transx3,transy3,pp3,pp3d,pp3dd,minx3,slope3r,ratio3r);
dist=@(xz,yz,point)distfunc(xz,yz,point,theta1,theta2,theta3,dist1,dist2,dist3);

xs=linspace(-.02,.2,numpts);
ys=linspace(-stepsize*numlayer,stepsize*numlayer,2*numlayer+1);
[xs,ys]=meshgrid(xs,ys);
[xs,ys]=rotatepts(xs(:),ys(:),(270-angle2)*pi/180,0,0);
phimap1=zeros(total,3);
points1=zeros(total,2);
points={};
phimap={};
for z=1:total
    x=xs(z);
    y=ys(z);
    points1(z,1)=x;
    points1(z,2)=y;
    [xr1,yr1]=rotatepts(x,y,-theta1,0,0);
    [xr2,yr2]=rotatepts(x,y,-theta2,0,0);
    [xr3,yr3]=rotatepts(x,y,theta3,0,0);

    dmin=distance(xr1,yr1,transx1,transy1,pp1,pp1d,pp1dd,minx1,slope1r,ratio1r);
    phimap1(z,1)=dmin;
    dmin=distance(xr2,yr2,transx2,transy2,pp2,pp2d,pp2dd,minx2,slope2r,ratio2r);
    phimap1(z,2)=dmin;
    dmin=distance(xr3,yr3,transx3,transy3,pp3,pp3d,pp3dd,minx3,slope3r,ratio3r);
    phimap1(z,3)=dmin;
end

points{1}=points1;
phimap{1}=phimap1;

if single
  
  num=0;
  numpts=(numlayer*2+1); 
  return
end

xs=linspace(-.02,.2,numpts);
ys=linspace(-stepsize*numlayer,stepsize*numlayer,2*numlayer+1);
[xs,ys]=meshgrid(xs,ys);
[xs,ys]=rotatepts(xs(:),ys(:),(angle3-90)*pi/180,0,0);
phimap2=zeros(total,3);
points2=zeros(total,2);
for z=1:total
    x=xs(z);
    y=ys(z);
    points2(z,1)=x;
    points2(z,2)=y;
    [xr1,yr1]=rotatepts(x,y,-theta1,0,0);
    [xr2,yr2]=rotatepts(x,y,-theta2,0,0);
    [xr3,yr3]=rotatepts(x,y,theta3,0,0);

    dmin=distance(xr1,yr1,transx1,transy1,pp1,pp1d,pp1dd,minx1,slope1r,ratio1r);
    phimap2(z,1)=dmin;
    dmin=distance(xr2,yr2,transx2,transy2,pp2,pp2d,pp2dd,minx2,slope2r,ratio2r);
    phimap2(z,2)=dmin;
    dmin=distance(xr3,yr3,transx3,transy3,pp3,pp3d,pp3dd,minx3,slope3r,ratio3r);
    phimap2(z,3)=dmin;

end

points{2}=points2;
phimap{2}=phimap2;

xs=linspace(-.02,.2,numpts);
ys=linspace(-stepsize*numlayer,stepsize*numlayer,2*numlayer+1);
[xs,ys]=meshgrid(xs,ys);
[xs,ys]=rotatepts(xs(:),ys(:),-pi/2,0,0);
phimap3=zeros((numlayer*2+1)*numpts,3);
points3=zeros((numlayer*2+1)*numpts,2);

for z=1:total
    x=xs(z);
    y=ys(z);
    points3(z,1)=x;
    points3(z,2)=y;
    [xr1,yr1]=rotatepts(x,y,-theta1,0,0);
    [xr2,yr2]=rotatepts(x,y,-theta2,0,0);
    [xr3,yr3]=rotatepts(x,y,theta3,0,0);

    dmin=distance(xr1,yr1,transx1,transy1,pp1,pp1d,pp1dd,minx1,slope1r,ratio1r);
    phimap3(z,1)=dmin;
    dmin=distance(xr2,yr2,transx2,transy2,pp2,pp2d,pp2dd,minx2,slope2r,ratio2r);
    phimap3(z,2)=dmin;
    dmin=distance(xr3,yr3,transx3,transy3,pp3,pp3d,pp3dd,minx3,slope3r,ratio3r);
    phimap3(z,3)=dmin;

end

points{3}=points3;
phimap{3}=phimap3;

num=ceil(.02/stepsize)+1;
xs=linspace(-.01,.01,num);
phimap4=zeros(num^2,3);
points4=zeros(num^2,2);
z=1;
for i=1:num
    for j=1:num
        x=xs(i);
        y=xs(j);
    
        points4(z,1)=x;
        points4(z,2)=y;
        [xr1,yr1]=rotatepts(x,y,-theta1,0,0);
        [xr2,yr2]=rotatepts(x,y,-theta2,0,0);
        [xr3,yr3]=rotatepts(x,y,theta3,0,0);
        
        dmin=distance(xr1,yr1,transx1,transy1,pp1,pp1d,pp1dd,minx1,slope1r,ratio1r);
        phimap4(z,1)=dmin;
        dmin=distance(xr2,yr2,transx2,transy2,pp2,pp2d,pp2dd,minx2,slope2r,ratio2r);
        phimap4(z,2)=dmin;
        dmin=distance(xr3,yr3,transx3,transy3,pp3,pp3d,pp3dd,minx3,slope3r,ratio3r);
        phimap4(z,3)=dmin;
        z=z+1;
    end
end

points{4}=points4;
phimap{4}=phimap4;
numpts=(numlayer*2+1);
end
function dmin=distance(xr,yr,transx,transy,pp,ppd,ppdd,minx,sloper,ratior)
    if xr<=minx
        dmin=abs(yr-xr*-sloper)*ratior;
        if yr>xr*-sloper
           dmin=dmin*-1; 
        end
    elseif xr>=-minx
        dmin=abs(yr-xr*sloper)*ratior;
        if yr>xr*sloper
           dmin=dmin*-1; 
        end
    else
        [x0]=argmin(xr,yr,transx,transy);
        [dmin,~]=distancetospline(xr,yr,x0,pp,ppd,ppdd);
    end

end

function [h,hx,hy,hxx,hxy,hyy]=distanceall(xr,yr,theta,transx,transy,pp,ppd,ppdd,minx,sloper,ratior)
    hyy=0;
    hxy=0;
    hxx=0;
    if xr<=minx
        h=(yr-xr*-sloper)*ratior;
        hx=(-cos(theta)*(-sloper)+sin(theta))/sqrt(1+(-sloper)^2);
        hy=(sin(theta)*(-sloper)+cos(theta))/sqrt(1+(-sloper)^2);
    elseif xr>=-minx
        h=(yr-xr*sloper)*ratior;
        hx=(-cos(theta)*(sloper)+sin(theta))/sqrt(1+(sloper)^2);
        hy=(sin(theta)*(sloper)+cos(theta))/sqrt(1+(sloper)^2);
    else
        [x0]=argmin(xr,yr,transx,transy);
        [h,hx,hy,hxx,hxy,hyy]=distancetosplineall(xr,yr,x0,theta,pp,ppd,ppdd);
    end

end

function [x0]=argmin(sx,sy,x,y)
    dist=(x-sx).^2+(y-sy).^2;
    [~,I]=min(dist);
    x0=x(I);
end