function val=convguass(x,y,px,py,DT,top)
n=size(x,2);
intergrand=(top*erf((y-py)/(2*sqrt(DT)))+1).*exp(-(x-px).^2/(4*DT))/(4*sqrt(pi)*sqrt(DT));
val=((x(2)-x(1))*intergrand(1)+sum((x(3:n)-x(1:n-2)).*intergrand(2:n-1))+(x(n)-x(n-1))*intergrand(n))/2;

