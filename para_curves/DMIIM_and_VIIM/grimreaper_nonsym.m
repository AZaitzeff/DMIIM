function [xl,yl,xr,yr] = grimreaper_nonsym(nt,dt,ren,angle1,angle2,gamma,mu,np,ruuthfunc,ruuthfunccenter,ratio)


ylim=.5;
S1=gamma(1);
S2=gamma(2);
S3=gamma(3);

mu1=mu(1);
mu2=mu(2);
mu3=mu(3);

[xl,yl,xr,yr,xb,yb,~]=VIIM_initialdata_nonsym(angle1,angle2,mu,np,ylim,0);
for t=1:nt
    
    x1=[xl xr(end-1:-1:1)];
    y1=[yl yr(end-1:-1:1)];
    
    x2=[xl xb(2:end)];
    y2=[yl yb(2:end)];
    
    x3=[xr xb(2:end)];
    y3=[yr yb(2:end)];



    [x1,y1] = kevolvefullyimplict(dt,ren,x1,y1,S1,mu1,0,0,0);
    [x2,y2] = kevolvefullyimplict(dt,ren,x2,y2,S2,mu2,4,0,0);
    [x3,y3] = kevolvefullyimplict(dt,ren,x3,y3,S3,mu3,4,0,0);
    if any(isnan([x1,y1,x2,y2,x3,y3])) || any(abs([x1,y1,x2,y2,x3,y3])>.5e2)
      save(['data/errorkev' num2str(angle1) num2str(np) num2str(round(mu1*100)) '.mat'],'x1','y1','x2','y2','x3','y3','S1','S2','S3',...
          'angle1','angle2','dt','ren','t','xl','yl','xr','yr','xb','yb')
      break
    end
    try
        [xl,yl,xr,yr,xb,yb]=ruuth_x_spline_unequal(x1,y1,x2,y2,x3,y3,xl,xr,xl(end),yl(end),angle1,angle2,ruuthfunc,ruuthfunccenter,ratio);
    catch
      save(['data/errorkev' num2str(angle1) num2str(np) num2str(round(mu1*100)) '.mat'],'x1','y1','x2','y2','x3','y3','S1','S2','S3',...
          'angle1','angle2','dt','ren','t','xl','yl','xr','yr','xb','yb')
      break
    end
    if mod(t,100)==0
       %DT=t*dt*ren;
       %save(strcat('data/interb',num2str(z)),'x','y','DT')
       %z=z+1;
       strcat(num2str(t/nt*100),' percent done')
   end
end
end
