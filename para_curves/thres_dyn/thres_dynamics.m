function [x,y] = thres_dynamics(angle,N,nt,DT,ST,boost)

[x,y,~]=VIIM_initialdata(angle,N,0,0);

for t=1:nt
    [x,y]=thres_dynamics_one_step(x,y,DT,ST,boost);
end
end