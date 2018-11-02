function [phi]=dictmapping2(phi,h,level,FLAG,map,funcs,DMIIM)
    

    [~,c]=size(phi);
    [n,m]=size(phi{1});
    all=1:c;
    redists=zeros(c,n,m);
    for i=1:c
       temp=reshape(phi{i}, [n m]);

       [temp,~,~] = redistz(temp,level,FLAG,h,h);% replace with own redistancing code
       redists(i,:,:)=temp;

    end
    maxvals=max(redists,[],1);
    
    for i=1:c 
        rest=setdiff(all,i);
        tempmaxvals=max(redists(rest,:,:),[],1);
        temp=redists(i,:,:)-tempmaxvals;
        phi{i}=reshape(temp,[m n]);
    end
    
    if DMIIM
    indices=find(maxvals<2.5*h);
    redistsflat=reshape(redists,[c n*m]);
    

    totalnum=numel(indices);
    for indp=1:totalnum
        [k,l]=ind2sub([n m],indices(indp));
        vals=redistsflat(:,indices(indp));
        [~,I]=sort(vals,'descend');
        key=sort(I(1:3))';
        op=map(mat2str(key));
        [lsv]=funcs{op.index}(vals(op.order));
        for i=1:3
            phi{op.order(i)}(k,l)=lsv(i);
        end
        val=min(lsv);
        for j=4:min(5,c)
            phi{I(j)}(k,l)=val;
        end
    end
    
    end


    for i=1:c

      temp=reshape(phi{i}, [n m]);
      [temp,~,~] = redistz(temp,level,FLAG,h,h);% replace with own redistancing code
      phi{i}=temp;

    end
end