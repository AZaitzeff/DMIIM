%Generates the data for tables 2 and 3
addpath('../utility')
addpath('../para_curves/DMIIM_and_VIIM')
oT1=18/512;

num=4;
for j = 1:2*num
    i=mod(j-1,num)+1;
    DT=1/(2^(12+i));
    ren=2^(12);
    dt=DT/ren;
    eps=1/2^(7+(i-1)/4);
    if j<=num
        np=2048;
        [~,~,~,~]=testanyangleFT(90,2,oT1,dt,ren,np,eps,0,0,strcat('para90',num2str(12+i),'DT',num2str(11),'np',num2str(27+i),'o4epsl'));
    else
        np=2048;
        [~,~,~,~]=testanyangleFT(120,2,oT1,dt,ren,np,eps,0,0,strcat('para120',num2str(12+i),'DT',num2str(11),'np',num2str(27+i),'o4epsl'));
    end
end