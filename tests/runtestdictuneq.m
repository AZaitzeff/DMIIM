%Generates the data for tables 12 and 13
addpath('../utility')
addpath('../para_curves/DMIIM_and_VIIM')


if exist('dict/dict75135Dtonesmu.mat', 'file')==0
    savedictionary2(75,135,[1,1,1],[],'ones',0)
end
if exist('dict/dict75135Dtinvmu.mat', 'file')==0
angle1=75;
angle2=135;
angle1r=angle1/180*pi;
angle2r=angle2/180*pi;
angle3r=2*pi-angle1r-angle2r;
mu=zeros(1,3);
mu(1)=1/(sin(angle3r)+sin(angle2r)-sin(angle1r));
mu(2)=1/(sin(angle3r)+sin(angle1r)-sin(angle2r));
mu(3)=1/(sin(angle1r)+sin(angle2r)-sin(angle3r));
savedictionary2(75,135,mu,[],'inv',0)
end

num=7;
for j = 1:2*num
    i=mod(j-1,num)+1;
    T=18/512;

    if j<=num
        np=ceil(2^(9+1/2+i/2));
        DT=1/(2^(9+i));
        ren=ceil(2^(9+1/2+i/2));
        dt=DT/ren;
        testanyangleFT_nonsym(75,135,0,1,DT,T,ren,np,0,strcat('dict75135',num2str(9+i),'DT',num2str(1),'np',num2str(1),'rencn'),'ones');
    elseif j<=2*num
        np=ceil(2^(9+1/2+i/2));
        DT=1/(2^(9+i));
        ren=ceil(2^(9+1/2+i/2));
        dt=DT/ren;
        testanyangleFT_nonsym(75,135,0,0,DT,T,ren,np,0,strcat('dict75135',num2str(9+i),'DT',num2str(1),'np',num2str(1),'reninv'),'inv');
    end
end