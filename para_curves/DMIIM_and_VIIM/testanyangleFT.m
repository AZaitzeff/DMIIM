function [x,y,truey,error]=testanyangleFT(angle,ngamma,oT,dt,ren,np,eps,mode,factor,name,impl,nameload)
    if nargin()<11
        impl=1;
        nameload='';
    end
        
    nbhrcor=.04;
    DT=ren*dt;
    %[gamma]=findgamFT(angle,DT,eps);
    %if numel(gam)>1
    %    gamma=gam;
    %end
    angler=angle/180*pi;
    gam=sin(angler)/sin((2*pi-angler)/2);
    gamma=[2-gam,gam,gam];
    %gamma=ngamma;
    nt=ceil(oT/(DT));
    T=nt*DT;
    tic;
    [x,y] = grimreaper_sym(nt,dt,ren,eps,angle,gamma,np,mode,impl,factor,nameload);
    toc;
    if mode==1 || mode==3
        [~,truey]=line_initialdata(angle,np);
    else
        [~,truey,~]=VIIM_initialdata(angle,np,eps,T);
    end
        %truey=real((1/coef)*log(cos(coef*x))-coef*T);
    difference=(truey-y);
    mask=abs((x-.25))>nbhrcor;
    x1=x(mask);
    dif=difference(mask);
    
    if mode==2
        trueareadif=abs(cot(angle*pi/180/2)*x1(end)^2);
    else
        trueareadif=(abs(truey(1))*(x1(end)));
    end
    
    error=trap_per(dif,x1)/trueareadif;
    
    save(strcat('data/',name,'.mat'),'x','y','truey','dt','nt','eps','ren','T','gamma','error')

end