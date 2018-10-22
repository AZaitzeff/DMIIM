function [error]=testanyangleFT_nonsym(angle1,angle2,gam,muoption,DT,T,ren,np,eps,name,loadname)



    if numel(gam)>1
        gamma=gam;
    else
        angle1r=angle1/180*pi;
        angle2r=angle2/180*pi;
        angle3r=2*pi-angle1r-angle2r;
        gamma=zeros(1,3);
        gamma(1)=(sin(angle3r)+sin(angle2r)-sin(angle1r));
        gamma(2)=(sin(angle3r)+sin(angle1r)-sin(angle2r));
        gamma(3)=(sin(angle1r)+sin(angle2r)-sin(angle3r));
    end
    
    if numel(muoption)>1
        mu=muoption;
    elseif muoption==0
        mu=zeros(1,3);
        mu(1)=1/gamma(1);
        mu(2)=1/gamma(2);
        mu(3)=1/gamma(3);
    else
        mu=[1,1,1];
    end
    %[~,~,~,~,~,~,vel]=VIIM_initialdata_nonsym(angle1,angle2,mu,10,.5,0);
    %T=dist/vel;
    nt=T/DT;
    dt=DT/ren;
    ratio=sqrt(2^(11)*DT);
    
    if strcmp(loadname,'')
        ruuthfunc=@(p,a)voronoireconst(p,a);
        ruuthfunccenter=@(p)voronoireconstcenter(p);
    else
        vars=load(['dict/dict' num2str(angle1) num2str(angle2) 'Dt' loadname 'mu.mat']);
        ruuthfunc=@(p,a)dictmapreconstz(p,a,vars.points,vars.phimap,vars.numpts,vars.slopes,ratio,vars.dist);
        ruuthfunccenter=@(p)dictmapreconstcenter(p,vars.points{4},vars.phimap{4},vars.num,ratio,vars.dist);
    end
    ti=cputime;
    [xl,yl,xr,yr] = grimreaper_nonsym(nt,dt,ren,angle1,angle2,gamma,mu,np,ruuthfunc,ruuthfunccenter,ratio);
    e=cputime-ti
    
    
    [error,trueyl,trueyr]=error_unequal(xl,yl,xr,yr,angle1,angle2,mu,T);
    
    save(strcat('data/',name,'.mat'),'xl','yl','trueyl','xr','yr','np','trueyr','angle1','angle2','dt','nt','eps','ren','DT','T','gamma','mu','error')

end