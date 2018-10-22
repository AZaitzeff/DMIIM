function [phi]=VIIMwDM2(phi,gamma,h,dt,eps,nt,ren,map,width,FLAG,funcs,DMIIM)
    
    for iter=1:nt
        phi=levelsetstepsign2(phi,gamma,h,dt,ren,eps);
        phi=dictmapping2(phi,h,width,FLAG,map,funcs,DMIIM);
    end
end

