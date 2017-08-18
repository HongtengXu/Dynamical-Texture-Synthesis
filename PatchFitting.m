function out=PatchFitting(overlap, patch, step, pr, L, mask)


over=2*pr+1-step;

m11=mask(1:over,1:over);
m12=mask(1:over,end-over+1:end);
m21=mask(end-over+1:end,1:over);
m22=mask(end-over+1:end,end-over+1:end);

LL=zeros(2);
LL(1,1)=double(sum(m11(:))==over^2);
LL(1,2)=double(sum(m12(:))==over^2);
LL(2,1)=double(sum(m21(:))==over^2);
LL(2,2)=double(sum(m22(:))==over^2);

if sum(LL(:))>1

    mleft=mask(:,1:over);
    mright=mask(:,end-over+1:end);
    mup=mask(1:over,:);
    mdown=mask(end-over+1:end,:);

    res=sum(abs(overlap-patch.*repmat(mask,[1,1,L])),3)+0.01;
    resleft=res(:,1:over);
    resright=res(:,end-over+1:end);
    resup=res(1:over,:);
    resdown=res(end-over+1:end,:);


    if sum(mleft(:))==( (2*pr+1)*over )
        C = mincut(resleft,0);
        mask(:,1:over)=mask(:,1:over) + double(C <= 0);
    end

    if sum(mright(:))==( (2*pr+1)*over )
        C = mincut(resright,0);
        mask(:,end-over+1:end)=mask(:,end-over+1:end) + double(C >= 0);
    end



    if sum(mup(:))==( (2*pr+1)*over )
        C = mincut(resup,1);
        mask(1:over,:)= mask(1:over,:) + double(C <= 0);
    end


    if sum(mdown(:))==( (2*pr+1)*over )
        C = mincut(resdown,1);
        mask(end-over+1:end,:)= mask(end-over+1:end,:) + double(C >= 0);
    end

else
    

    res=sum(abs(overlap-patch.*repmat(mask,[1,1,L])),3)+0.01;
    resleft=res(1:over,1:over);
    resright=res(1:over,end-over+1:end);
    resup=res(end-over+1:end,1:over);
    resdown=res(end-over+1:end,end-over+1:end);


    if LL(1,1)==1
        C1 = mincut(resleft,0);
        C2 = mincut(resleft,1);
        mask(1:over,1:over)=mask(1:over,1:over) + double(C1 <= 0) + double(C2 <= 0);
    end


    if LL(1,2)==1
        C1 = mincut(resright,0);
        C2 = mincut(resright,1);
        mask(1:over,end-over+1:end)= mask(1:over,end-over+1:end) + double(C1 >= 0) + double(C2 <= 0);
    end



    if LL(2,1)==1
        C1 = mincut(resup,0);
        C2 = mincut(resup,1);
        mask(end-over+1:end,1:over)= mask(end-over+1:end,1:over) + double(C1 <= 0) + double(C2 >= 0);
    end


    if LL(2,2)==1
        C1 = mincut(resdown,0);
        C2 = mincut(resdown,1);
        mask(end-over+1:end,end-over+1:end)=mask(end-over+1:end,end-over+1:end) + double(C1 >= 0) + double(C2 >= 0);
    end
    
end

mask=double(mask>1);
out=overlap.*repmat(mask,[1,1,L])+patch.*repmat(1-mask,[1,1,L]);