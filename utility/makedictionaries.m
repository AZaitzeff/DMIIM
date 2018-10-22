for i=3:4
if i==1
	savedictionary2(90,135,[1,1,1],[2-sqrt(2),sqrt(2),sqrt(2)],'VIIM',1)
elseif i==2
savedictionary2(120,120,[1,1,1],[1,1,1],'normal',1)
elseif i==3
savedictionary2(75,135,[1,1,1],[],'ones',0)
elseif i==4
angle1=75;
angle2=135;
angle1r=angle1/180*pi;
angle2r=angle2/180*pi;
angle3r=2*pi-angle1r-angle2r;
mu=zeros(1,3);
mu(1)=1/(sin(angle3r)+sin(angle2r)-sin(angle1r));
mu(2)=1/(sin(angle3r)+sin(angle1r)-sin(angle2r));
mu(3)=1/(sin(angle1r)+sin(angle2r)-sin(angle3r));
savedictionary2(75,135,mu,[],'inv',0)
elseif i==6
savedictionarygrid(120,120,[1,1,1],'ones')
elseif i==7
savedictionarygrid(90,135,[1,1,1],'ones')
elseif i==8
savedictionarygrid(120,120,[2-sqrt(2),2-sqrt(2),sqrt(2)],'spec')
end
end