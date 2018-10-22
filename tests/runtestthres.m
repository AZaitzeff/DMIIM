%Generates the data for tables 6 and 7
addpath('../utility')
addpath('../para_curves/thres_dyn')

oT1=18/512;
num=7;
parfor j = 1:2*num
    i=mod(j-1,num)+1;
    if j<=num
        DT=1/(2^(10+i));
        np=2^(11+i);
        nt=round(oT1/DT);
        [~,~,~,~]=testanyanglethes(90,oT1,nt,np,0,strcat('para90',num2str(10+i),'DT',num2str(11+i),'npthres'));
    else
        DT=1/(2^(11+i));
        np=2^(11+i);
        nt=round(oT1/DT);
        [~,~,~,~]=testanyanglethes(90,oT1,nt,np,0,strcat('para90',num2str(11+i),'DT',num2str(11+i),'npthres'));
    end
end