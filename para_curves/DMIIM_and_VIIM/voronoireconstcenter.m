function [lsv]=voronoireconstcenter(point)
    lsv=(point(1)-point(2))^2+(point(2)-point(3))^2+(point(1)-point(3))^2;
end