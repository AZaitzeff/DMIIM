N=512;
nt=80;
%N=128;
%nt=10;

width=N;   
addpath('../utility/redist/')
addpath('../utility/')
addpath('../level_set_DMIIM/')
if exist('dict/dict120grgrid.mat', 'file')==0
    savedictionarygrid(120,120,[sqrt(2),2-sqrt(2),2-sqrt(2)],'120gr')
end
if exist('dict/dict90grgrid.mat', 'file')==0
    savedictionarygrid(90,135,[2-sqrt(2),sqrt(2),sqrt(2)],'90gr')
end

angle1=90;
angle2=135;
angle1r=angle1/180*pi;
angle2r=angle2/180*pi;
angle3r=2*pi-angle1r-angle2r;
%gamma=zeros(1,4);
%gamma(1)=(sin(angle3r)+sin(angle2r)-sin(angle1r));
%gamma(2)=(sin(angle3r)+sin(angle1r)-sin(angle2r));
%gamma(3)=(sin(angle1r)+sin(angle2r)-sin(angle3r));
%gamma(4)=gamma(1);
gamma(1)=2-sqrt(2);
gamma(2)=sqrt(2);
gamma(3)=sqrt(2);
gamma(4)=2-sqrt(2);

%N=128;
%nt=20;
[u,v,w,vel] = grimreaper(angle1,angle2,0,N,[1,1,1]);
h= 1/N;

z=u;
u(end,:)=-abs(u(end,:));
u(N/2+1:end-1,:)=-1;
z(1,:)=-abs(z(1,:));
z(2:N/2,:)=-1;
phi={};
[temp,~,~] = redistz(u,width,1,h,h);% replace with own redistancing code
phi{1}=temp;
[temp,~,~] = redistz(v,width,1,h,h);% replace with own redistancing code
phi{2}=temp;
[temp,~,~] = redistz(w,width,1,h,h);% replace with own redistancing code
phi{3}=temp;
[temp,~,~] = redistz(z,width,1,h,h);% replace with own redistancing code
phi{4}=temp;
T=.11/vel;
dt=h^2/(5*max(gamma));
DT=T/nt;
ren=ceil(DT/dt);
dt=DT/ren;

map = containers.Map;
key=[1 2 3];
op.index=1;
op.order=key;
map(mat2str(key))=op;
key=[2 3 4];
op.index=1;
op.order=[4 3 2];
map(mat2str(key))=op;
key=[1 3 4];
op.index=2;
op.order=[3 4 1];
map(mat2str(key))=op;
key=[1 2 4];
op.index=2;
op.order=[2 4 1];
map(mat2str(key))=op;

funcs={};
angles=[90,120];
for i=1:2
    vars=load(['dict/dict' num2str(angles(i)) 'grgrid.mat']);
    ratio=sqrt(2^(11)*DT);
    funcs{i}= @(p) dictmaprecgrid(p,vars.points,vars.phimap,vars.num,...
        vars.numpts*vars.num,vars.slopes,ratio,vars.dist);
end
%%
ti=cputime;

nind={};
nind{1}=phi;
for j=1:nt
    [phi]=VIIMwDM2(phi,gamma,h,dt,0,1,ren,map,width,1,funcs,1);
    nind{j+1}=phi;
end
error=0;

%error=0;
save(strcat('data/','grim','.mat'),'angle1','angle2','N','h','dt','nt','ren','T','vel','gamma','nind')
e=cputime-ti