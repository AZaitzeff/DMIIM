%Generates data for tables 1o and 11

addpath('../utility')
addpath('../para_curves/DMIIM_and_VIIM')

oT1=18/512;
num=7;
if exist('dict/dict90135DtVIIMmu.mat', 'file')==0
    savedictionary2(90,135,[1,1,1],[2-sqrt(2),sqrt(2),sqrt(2)],'VIIM',1)
end
if exist('dict/dict90135Dtnormalmu.mat', 'file')==0
    savedictionary2(120,120,[1,1,1],[1,1,1],'normal',1)
end



for j = 1:2*num
    i=mod(j-1,num)+1;
    oT1=18/512;
    if j<=num
        np=ceil(2^(9+1/2+i/2));
        DT=1/(2^(9+i));
        ren=ceil(2^(9+1/2+i/2));
        dt=DT/ren;
        [~,~,~,~]=testanyangleFT(90,2,oT1,dt,ren,np,0,6,0,strcat('dict90',num2str(9+i),'DT',num2str(1),'np',num2str(1),'ren'),1,'VIIM');
    elseif j<=2*num
        np=ceil(2^(9+1/2+i/2));
        DT=1/(2^(9+i));
        ren=ceil(2^(9+1/2+i/2));
        dt=DT/ren;
        [~,~,~,~]=testanyangleFT(120,2,oT1,dt,ren,np,0,6,0,strcat('dict120',num2str(9+i),'DT',num2str(1),'np',num2str(1),'ren'),1,'normal');
    end
end
