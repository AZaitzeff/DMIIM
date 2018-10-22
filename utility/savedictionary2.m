function savedictionary2(angle1,angle2,mu,gamma,name,single)
    angle1r=angle1/180*pi;
	angle2r=angle2/180*pi;
	angle3r=2*pi-angle1r-angle2r;
    if numel(gamma)<2
        gamma=zeros(1,3);
        gamma(1)=(sin(angle3r)+sin(angle2r)-sin(angle1r))*mu(1);
        gamma(2)=(sin(angle3r)+sin(angle1r)-sin(angle2r))*mu(2);
        gamma(3)=(sin(angle1r)+sin(angle2r)-sin(angle3r))*mu(3);
    end
    [points,phimap,slopes,numpts,num,dist]=makedictionary2(angle1,angle2,2^-11,gamma,5e-4,21,single);
    save(['dict/dict' num2str(angle1) num2str(angle2) 'Dt' name 'mu.mat'],'points','phimap','slopes','numpts','num','dist');
end