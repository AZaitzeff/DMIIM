%Generates the data for table 4 and 5
addpath('../utility')
addpath('../para_curves/DMIIM_and_VIIM')


oT1=18/512;
num=4;
for j = 1:4*num
    i=mod(j-1,num)+1;
    DT=1/(2^(12+i));
    ren=2^(12);
    dt=DT/ren;
    if j<=num
        eps=2*sqrt(DT);
        np=2048;
        [~,~,~,~]=testanyangleFT(90,2,oT1,dt,ren,np,eps,0,0,strcat('para90',num2str(12+i),'DT',num2str(11),'np',num2str(2),'epsh'));
    elseif j<=2*num
        eps=2*sqrt(DT);
        np=2048;
        [~,~,~,~]=testanyangleFT(120,2,oT1,dt,ren,np,eps,0,0,strcat('para120',num2str(12+i),'DT',num2str(11),'np',num2str(2),'epsh'));
    elseif j<=3*num
        eps=4*sqrt(DT);
        np=2048;
        [~,~,~,~]=testanyangleFT(90,2,oT1,dt,ren,np,eps,0,0,strcat('para90',num2str(12+i),'DT',num2str(11),'np',num2str(4),'epsh'));
    else
        eps=4*sqrt(DT);
        np=2048;
        [~,~,~,~]=testanyangleFT(120,2,oT1,dt,ren,np,eps,0,0,strcat('para120',num2str(12+i),'DT',num2str(11),'np',num2str(4),'epsh'));
    end
end
