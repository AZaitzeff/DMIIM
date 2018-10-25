rng(104)
addpath('../utility/redist/')
addpath('../utility/')
addpath('../level_set_DMIIM/')
if exist('dict/dict120grid.mat', 'file')==0
    savedictionarygrid(120,120,[1,1,1],'120')
end
if exist('dict/dict1202grid.mat', 'file')==0
    savedictionarygrid(120,120,[sqrt(2)-1,sqrt(2)-1,sqrt(2)-1],'1202')
end
if exist('dict/dict90grid.mat', 'file')==0
savedictionarygrid(90,135,[sqrt(2)-1,1,1],'90')
end
if exist('dict/dict146grid.mat', 'file')==0
ang=2*acos((sqrt(2)-1)/sqrt(2))*180/pi;
savedictionarygrid(ang,180-ang/2,[1,sqrt(2)-1,sqrt(2)-1],'146')
end

DMIIM=1;
flag=1;

map = containers.Map;
jump=3;
cases=3;
num=4;
total=num^2;
for i=1:total-2
    for j=(i+1):total-1
        for k=(j+1):total
            key=[i j k];
            st=[(mod(i,jump)==0) (mod(j,jump)==0) (mod(k,jump)==0)];
            totalv=sum(st);
            op.index=totalv+1;
            op.order=key;
            if totalv==1
                ind=find(st==1);
                temp=op.order(ind);
                op.order(ind)=op.order(1);
                op.order(1)=temp;
            elseif totalv==2
                ind=find(st==0);
                temp=op.order(ind);
                op.order(ind)=op.order(1);
                op.order(1)=temp;
                
            end
            map(mat2str(key))=op;
        end
    end
end


N=200;
width=N;
x=linspace(0,1-1/N,N);
h=x(2)-x(1);
dt=h^2/10;
ren=100;
DT=dt*ren;
[X,Y]=meshgrid(x);
phi={};
points=zeros(total,2);
[xpts,ypts]=meshgrid(((1:num)-1/2)/num);
points(:,1)=xpts(:);
points(:,2)=ypts(:);
points=points+rand(total,2)*1/(2*num);
for z=1:total
    phi{z}=-(min(min((points(z,1)+1-X).^2,(points(z,1)-X).^2),(points(z,1)-1-X).^2)+...
        min(min((points(z,2)+1-Y).^2,(points(z,2)-Y).^2),(points(z,2)-1-Y).^2));
end
phitemp=phi;
for z=1:total
    maxmatrix=ones(N,N)*-inf;
    for i=1:total
        if i~=z
            maxmatrix=max(maxmatrix,phitemp{i});
        end
        
    end
    [temp,~,~] = redistz(phitemp{z}-maxmatrix,width,flag,1/N,1/N);% replace with own redistancing code
    phi{z}=temp;
    
end

angles=[120,90,146,120];
funcs={};
for i=1:4
    if i<4
    vars=load(['../poly/ruuthdict/dict' num2str(angles(i)) 'mgrid.mat']);
    else
        vars=load(['../poly/ruuthdict/dict' num2str(angles(i)) 'm2grid.mat']);
    end
    ratio=sqrt(2^(11)*DT);
    funcs{i}= @(p) dictmaprecgrid(p,vars.points,vars.phimap,vars.num,...
        vars.numpts*vars.num,vars.slopes,ratio,vars.dist);
end
%%
gamma=zeros(1,total);
for i=1:total
    if i==3
        gamma(i)=sqrt(2)-1;
    else
        gamma(i)=1;
    end
end
phis={};
phis{1}=phi;
for t=1:100
    [phi]=VIIMwDM2(phi,gamma,h,dt,0,1,ren,map,width,flag,funcs,DMIIM);
    phis{t+1}=phi;
    t
    save(['data/2dsqr' num2str(flag) 'DMIIM' num2str(DMIIM) '.mat'],'phis','dt','ren');
end