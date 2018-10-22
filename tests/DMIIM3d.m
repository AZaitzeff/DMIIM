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
width=100;
map = containers.Map;
jump=3;
cases=3;
for i=1:6
    for j=(i+1):7
        for k=(j+1):8
            key=[i j k];
            st=[(mod(i,jump)==0) (mod(j,jump)==0) (mod(k,jump)==0)];
            total=sum(st);
            op.index=mod(total,cases)+1;
            op.order=key;
            if total==1
                ind=find(st==1);
                temp=op.order(ind);
                op.order(ind)=op.order(1);
                op.order(1)=temp;
            elseif total==2
                ind=find(st==0);
                temp=op.order(ind);
                op.order(ind)=op.order(1);
                op.order(1)=temp;
                
            end
            map(mat2str(key))=op;
        end
    end
end

N=100;
x=linspace(0,1-1/N,N);
h=x(2)-x(1);
dt=h^2/10;
ren=100;
DT=dt*ren;
[X,Y,Z]=meshgrid(x);
num=2;
total=num^3;
points=zeros(total,3);
[xpts,ypts,zpts]=meshgrid(((1:num)-1/2)/num);
points(:,1)=xpts(:);
points(:,2)=ypts(:);
points(:,3)=zpts(:);
points=points+rand(total,3)*1/(2*num);
phi={};
z=1;

for z=1:total
    phi{z}=-(min(min((points(z,1)+1-X).^2,(points(z,1)-X).^2),(points(z,1)-1-X).^2)+...
        min(min((points(z,2)+1-Y).^2,(points(z,2)-Y).^2),(points(z,2)-1-Y).^2)+...
        min(min((points(z,3)+1-Z).^2,(points(z,3)-Z).^2),(points(z,3)-1-Z).^2));
end
phitemp=phi;
for z=1:total
    maxmatrix=ones(N,N,N)*-inf;
    for i=1:total
        if i~=z
            maxmatrix=max(maxmatrix,phitemp{i});
        end
        
    end
    temp=phitemp{z}-maxmatrix;
    %[temp,~,~] = redistz3(phitemp{z}-maxmatrix,width,flag,N,1/N,1/N,1/N);% replace next two lines with own redistancing code
    %temp=reshape(temp,[N N N]);
    phi{z}=temp*N*h;
    
end
%%
angles=[120,90,146];
funcs={};
for i=1:3
    vars=load(['dict' num2str(angles(i)) 'grid.mat']);
    ratio=sqrt(2^(11)*DT);
    funcs{i}= @(p) dictmaprecgrid(p,vars.points,vars.phimap,vars.num,...
        vars.numpts*vars.num,vars.slopes,ratio,vars.dist);
end
gamma=[1,1,sqrt(2)-1,1,1,sqrt(2)-1,1,1];
phis={};
phis{1}=phi;
for t=1:20
    [phi,~]=VIIMwDM3(phi,gamma,h,dt,1,ren,map,width,flag,funcs,DMIIM);
    phis{t+1}=phi;
    t
    save(['data/3dsqr' num2str(flag) 'DMMIM' num2str(DMIIM) '.mat'],'phis','dt','ren','-v7.3');
end

