
%Generates the data for table 6 and 7
addpath('../utility')
addpath('../para_curves/DMIIM_and_VIIM')

oT1=18/512;
num=7;
for j = 1:2*num
    i=mod(j-1,num)+1;
    DT=1/(2^(12+i));
    ren=2^(12);
    dt=DT/ren;
    if j<=num
        np=2048;
        [~,~,~,~]=testanyangleFT(90,2,oT1,dt,ren,np,0,0,0,strcat('para90',num2str(12+i),'DT',num2str(11),'np'));
    else
        np=2048;
        [~,~,~,~]=testanyangleFT(120,2,oT1,dt,ren,np,0,0,0,strcat('para120',num2str(12+i),'DT',num2str(11),'np'));
    end
end