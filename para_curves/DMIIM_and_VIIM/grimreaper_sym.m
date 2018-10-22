function [x,y] = grimreaper_sym(nt,dt,ren,eps,angle,gamma,N,mode,impl,factor,name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grim reaper test.

%         3

%   ---- \ / -----

%     1   |    2

%         | 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
FLAG=0;
tol=1e-6;
% Surface tensions and mobilities:

S1 = gamma(1); S2 = gamma(2);
ST=ones(1,3);

Dt=dt*ren;
flip=.5;
if mode==6
    
    
    ratio=sqrt(2^(11)*Dt);
    angle2=(360-angle)/2;
    vars=load(['dict/dict' num2str(angle) num2str(angle2) 'Dt' name 'mu.mat']);
    ruuthfunc=@(p,a)dictmapreconst(p,a,vars.points,vars.phimap,vars.numpts,vars.slopes,ratio,vars.dist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,y,xe]=VIIM_initialdata(angle,N,eps,0);

%z=1;
for t=1:nt
   if eps>2*tol
       [nx1,ny1]=epsilon_x_spline(x,y,xe,eps,-1);
       [nx2,ny2]=epsilon_x_spline(x,y,x,eps,1);
   else
        nx1=x;
        ny1=y;
        nx2=x;
        ny2=y;
   end
  
  
  M=floor((flip+ny1(end))/sqrt((ny1(end)-ny1(end-1))^2+(nx1(end)-nx1(end-1))^2));
  x2 = (0.25-eps)*ones(1,M+1);
  
  y2 =linspace(ny1(end),-flip,M+1);
  inx=zeros(1,N+M);
  iny=zeros(1,N+M);
  
  
  
  inx(1:N)=nx1;
  iny(1:N)=ny1;
  inx(N+1:N+M)=x2(2:M+1);
  iny(N+1:N+M)=y2(2:M+1);



  

  dir=0;%boundary condition
  
  
  if impl
      %[inx,iny] = kevolvecn(dt,ren,inx,iny,S1,1,4,dir,1/(N^2)*factor);
      %[nx2,ny2] = kevolvecn(dt,ren,nx2,ny2,S2,1,0,dir,1/(N^2)*factor);
      [inx,iny] = kevolvefullyimplict(dt,ren,inx,iny,S2,1,4,dir,1/(N^2)*factor);
      [nx2,ny2] = kevolvefullyimplict(dt,ren,nx2,ny2,S1,1,0,dir,1/(N^2)*factor);
  else
      [inx,iny] = kevolve(dt,ren,inx,iny,S2,1,4,dir,2);
      [nx2,ny2] = kevolve(dt,ren,nx2,ny2,S1,1,0,dir,2);
  end
  if any(isnan([inx,iny,nx2,ny2])) || any(abs([inx,iny,nx2,ny2])>.5e2)
      save(['data/errorkev' num2str(angle) num2str(N) '.mat'],'inx','iny','nx2','ny2','x','y','gamma','ST','angle','Dt','t')
      break
  end
  if mode==6
    [x,y]=ruuth_x_spline_nm(nx2,ny2,inx,iny,x,.25,y(end),angle,ruuthfunc,ratio);
    %[x,y,FLAG]=ruuth_x_spline(inx,iny,nx2,ny2,x,y(end),angle,points,phimap,xvalid);
  else
    [x,y,FLAG]=voronoi_x_spline(inx,iny,nx2,ny2,x,y(end),angle,ST);
  end

  
  if FLAG
      save(['data/errorvor' num2str(angle) num2str(N) '.mat'],'inx','iny','nx2','ny2','x','gamma','ST','angle','Dt','t')
      break
  end


   if mod(t,100)==0
       %DT=t*dt*ren;
       %save(strcat('data/interb',num2str(z)),'x','y','DT')
       %z=z+1;
       strcat(num2str(t/nt*100),' percent done')
   end
    

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%