function [phi]=dictmapping3(phi,h,map,width,flag,funcs,DMIIM)
    
    [~,c]=size(phi);
    [n,~,~]=size(phi{1});
    redists=zeros(c,n,n,n);
    all=1:c;
    for i=1:c
       temp=phi{i};
       temp(abs(temp)<1e-6)=0;
       [temp,~,~] = redistz3(temp,width,flag,n,1/n,1/n,1/n);% replace next two lines with own redistancing code
       temp=reshape(temp,[n n n]);
       redists(i,:,:,:)=temp*(h*n);
    end
    maxvals=max(redists,[],1);
    for i=1:c 
        rest=setdiff(all,i);
        tempmaxvals=max(redists(rest,:,:,:),[],1);
        temp=redists(i,:,:,:)-tempmaxvals;
        phi{i}=reshape(temp,[n n n]);
    end
    
    if flag==1
        indices=find(maxvals<2*h);
    else
        indices=find(maxvals<100*h);
    end
    redistsflat=reshape(redists,[c n^3]);
    totalnum=numel(indices);

    for indp=1:totalnum
        [j,k,l]=ind2sub([n n n],indices(indp));
        vals=redistsflat(:,indices(indp));
        [~,I]=sort(vals,'descend');
        key=sort(I(1:3))';
        op=map(mat2str(key));
        [lsv]=funcs{op.index}(vals(op.order));
        for i=1:3
            phi{op.order(i)}(j,k,l)=lsv(i);
        end
        val=min(lsv);
        for i=4:5
            phi{I(i)}(j,k,l)=val;
        end
    end
    %save('data/fir1.mat','phi')
    for i=1:c
       temp=phi{i};
       temp(abs(temp)<1e-6)=0;
       [temp,~,~] = redistz3(temp,width,flag,n,1/n,1/n,1/n);% replace next two lines with own redistancing code
       temp=reshape(temp,[n n n]);
       phi{i}=temp*(h*n);
    end
    %save('data/fir2.mat','phi')
end