function [x,y,truey,error]=testanyanglethes(angle,T,nt,np,boost,name)
    %do not use boost
    nbhrcor=.04;
    DT=T/nt;
    angler=angle/180*pi;
    gam=sin(angler)/sin((2*pi-angler)/2);
    ST=[[0,gam,1];[gam,0,1];[1,1,0]];
    T=nt*DT;
    tic;
    [x,y] = thres_dynamics(angle,np,nt,DT,ST,boost);
    toc;

    [~,truey,~]=VIIM_initialdata(angle,np,0,T);

    difference=(truey-y);
    mask=abs((x-.25))>nbhrcor;
    x1=x(mask);
    dif=difference(mask);
    

    trueareadif=(abs(truey(1))*(x1(end)));

    
    error=trap_per(dif,x1)/trueareadif;
    inferror=infnorm(x,y,T,angle);
    
    save(strcat('data/',name,'.mat'),'x','y','truey','DT','nt','T','ST','error','inferror')

end