function savedictionarygrid(angle1,angle2,mugamma,name)  
    [points,phimap,slopes,numpts,num,dist]=makedictionarygrid(angle1,angle2,2^-11,mugamma,1e-3,250);
    save(['dict/dict' name 'grid.mat'],'points','phimap','slopes','numpts','num','dist');
end