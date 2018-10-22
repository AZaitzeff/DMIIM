function [phi,phils]=VIIMwDM3(phi,gammamu,h,dt,nt,ren,map,width,flag,funcs,DMIIM)  
    addpath('../fast-redist/')
    addpath('../utility/')
    
    for iter=1:nt
        phils=levelsetstepsign3(phi,gammamu,h,dt,ren,width,flag);
        phi=dictmapping3(phils,h,map,width,flag,funcs,DMIIM);

    end
end
