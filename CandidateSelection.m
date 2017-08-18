function candidates=CandidateSelection(S1,ST,Samples,overlap,mask,K,T,flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% candidates selection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tmp1=Samples;
TmpT=Samples;
candidates=cell(T-2,1);
S=floor((T-2)/2);
if flag==1 % overlap is useless
    for t=1:S
        if isempty(Tmp1)
            Tmp1=Samples;
        end
        if isempty(TmpT)
            TmpT=Samples;
        end
        if t==1
            if size(Tmp1,2)<K
                Tmp1=Samples;
            end
            res1=sum((repmat(S1,[1,size(Tmp1,2)])-Tmp1).^2);
            [~,index] = sort(res1);           
            neighborhood = index(1:K);
            candidates{t}=Tmp1(:,neighborhood);
            Tmp1=Tmp1(:,index(K+1:end));

            if size(TmpT,2)<K
                TmpT=Samples;
            end
            resT=sum((repmat(ST,[1,size(TmpT,2)])-TmpT).^2);
            [~,index] = sort(resT);
            neighborhood = index(1:K);
            candidates{T-1-t}=TmpT(:,neighborhood);
            TmpT=TmpT(:,index(K+1:end));

        else

            cand1=candidates{t-1};
            candT=candidates{T-t};
            
            tmpcand1=[];       
            for i=1:size(cand1,2)
                if size(Tmp1,2)<size(cand1,2)
                    Tmp1=Samples;
                end
                res1=sum((repmat(cand1(:,i),[1,size(Tmp1,2)])-Tmp1).^2);
                resT=sum((repmat(candT(:,i),[1,size(Tmp1,2)])-Tmp1).^2);
                res = ( (T-1-(t-1))*res1+(t-1)*resT )./(T-1); 
                [~,index] = sort(res);
                tmpcand1=[tmpcand1,Tmp1(:,index(1))];
                Tmp1=Tmp1(:,index(2:end));
            end
            candidates{t}=tmpcand1;

            
            tmpcandT=[];
            for i=1:size(candT,2)
                if size(TmpT,2)<size(candT,2)
                    TmpT=Samples;
                end
                res1=sum((repmat(cand1(:,i),[1,size(TmpT,2)])-TmpT).^2);
                resT=sum((repmat(candT(:,i),[1,size(TmpT,2)])-TmpT).^2);
                res = ( (T-1-(t-1))*resT+(t-1)*res1 )./(T-1); 
                [~,index] = sort(res);
                tmpcandT=[tmpcandT,TmpT(:,index(1))];
                TmpT=TmpT(:,index(2:end));
            end
            candidates{T-1-t}=tmpcandT;
        end
    end
else
    
    X1=Samples.*repmat(mask,[1,size(Samples,2)]);
    XT=Samples.*repmat(mask,[1,size(Samples,2)]);
    for t=1:S
        if isempty(Tmp1) || ( size(Tmp1,2)<K)
            Tmp1=Samples;
            X1=Samples;
        end
        if isempty(TmpT) || (size(TmpT,2)<K)
            TmpT=Samples;
            XT=Samples;
        end
        
        if t==1
            resO=sum((repmat(overlap(:,t),[1,size(X1,2)])-X1).^2);
            res1=sum((repmat(S1,[1,size(Tmp1,2)])-Tmp1).^2);
            res=res1+size(Tmp1,1)/sum(mask(:))*resO;
            [~,index] = sort(res);
            Ktmp=min([K,size(Tmp1,2)]);
            neighborhood = index(1:Ktmp);
            candidates{t}=Tmp1(:,neighborhood);
            Tmp1=Tmp1(:,index(Ktmp+1:end));
            X1=X1(:,index(Ktmp+1:end));

            resO=sum((repmat(overlap(:,T-1-t),[1,size(XT,2)])-XT).^2);
            resT=sum((repmat(ST,[1,size(TmpT,2)])-TmpT).^2);
            res=resT+size(TmpT,1)/sum(mask(:))*resO;
            [~,index] = sort(res);
            Ktmp=min([K,size(TmpT,2)]);
            neighborhood = index(1:Ktmp);
            candidates{T-1-t}=TmpT(:,neighborhood);
            TmpT=TmpT(:,index(Ktmp+1:end));
            XT=XT(:,index(Ktmp+1:end));
            
        else
            
            cand1=candidates{t-1};
            candT=candidates{T-t};
            
            tmpcand1=[];       
            for i=1:size(cand1,2)
                if size(Tmp1,2)<size(cand1,2)
                    Tmp1=Samples;
                    X1=Samples;
                end
                resO=sum((repmat(overlap(:,t),[1,size(X1,2)])-X1).^2);
                res1=sum((repmat(cand1(:,i),[1,size(Tmp1,2)])-Tmp1).^2);
                resT=sum((repmat(candT(:,i),[1,size(Tmp1,2)])-Tmp1).^2);
                res = ( (T-1-(t-1))*res1+(t-1)*resT )./(T-1);               
                res=res+size(Tmp1,1)/sum(mask(:))*resO;
                [~,index] = sort(res);
                tmpcand1=[tmpcand1,Tmp1(:,index(1))];
                Tmp1=Tmp1(:,index(2:end));
                X1=X1(:,index(2:end));
            end
            candidates{t}=tmpcand1;

            tmpcandT=[];
            for i=1:size(candT,2)
                if size(TmpT,2)<size(candT,2)
                    TmpT=Samples;
                    XT=Samples;
                end
                resO=sum((repmat(overlap(:,T-1-t),[1,size(XT,2)])-XT).^2);
                res1=sum((repmat(cand1(:,i),[1,size(TmpT,2)])-TmpT).^2);
                resT=sum((repmat(candT(:,i),[1,size(TmpT,2)])-TmpT).^2);
                res = ( (T-1-(t-1))*resT+(t-1)*res1 )./(T-1);  
                res=res+size(TmpT,1)/sum(mask(:))*resO;
                [~,index] = sort(res);
                tmpcandT=[tmpcandT,TmpT(:,index(1))];
                TmpT=TmpT(:,index(2:end));
                XT=XT(:,index(2:end));
            end
            candidates{T-1-t}=tmpcandT;
            
        end
                       
    end
end

tmp=cell(T,1);
tmp{1}=S1;
tmp{T}=ST;
for t=2:T-1
    tmp{t}=candidates{t-1};
end
candidates=tmp;