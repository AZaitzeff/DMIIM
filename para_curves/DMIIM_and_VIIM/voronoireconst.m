function [lsv,phase]=voronoireconst(point,access)

if access==1
    lsv=point(1)-point(2);
    if lsv>=0
        phase=2;
    else
        phase=1;
    end
elseif access==2
    lsv=point(1)-point(3);
    if lsv>=0
        phase=3;
    else
        phase=1;
    end
else
    lsv=point(2)-point(3);
    if lsv>=0
        phase=3;
    else
        phase=2;
    end
end
lsv=abs(lsv);

end