function [W,X,Mask]=GraphCreatation(Candidates, sigmaN, K, d, Epsilon,method,beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gong D, Zhao X, Medioni G. Robust multiple manifolds structure learning
% [J]. arXiv preprint arXiv:1206.4624, 2012.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

% STEP1: COMPUTE PAIRWISE DISTANCES & FIND NEIGHBORS 
T=size(Candidates,1);
X=Candidates{1};
Mask=[];
for i=2:T
    X=[X,Candidates{i}];
    K1=size(Candidates{i-1},2);
    K2=size(Candidates{i},2);
    tmp=ones(K1,K2);
    Mask=blkdiag(Mask,tmp);
end
Mask=[zeros(size(Mask,1),1),Mask];
Mask=[Mask;zeros(1,size(Mask,2))];
Mask=Mask+Mask';
[~,N]=size(X);


X2 = sum(X.^2,1);
distance = repmat(X2,N,1)+repmat(X2',1,N)-2*(X'*X);
distance=distance.*Mask;
[~,index] = sort(distance);
neighborhood = index(2:(1+K),:);

if strcmp(method,'Eu')
    W=zeros(N);
    for i=1:N
        ind=find(distance(i,:)>0);
        for j=ind
            W(i,j)=distance(i,j)+Epsilon;
        end
    end
    if isnan(W(i,j))
        W(i,j)=0;
    end
else
    if strcmp(method,'Ag') 
        J=cell(N,1);
        for ii=1:N
            Xhat = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
            C = Xhat'*Xhat;                                        % local covariance
            S = diag(C);
            S = diag( 1./(1+sigmaN.*S)  ) ;
            [J{ii},~] = eigs(Xhat*(S*S')*Xhat',d);
        end
        W=zeros(N);
        for i=1:N
            ind=find(distance(i,:)>0);
            for j=ind
                [~,Angle]=prinAngles(J{i},J{j},'fast');
                Ang=sum(abs(Angle).^2);               
                W(i,j)=1/(Ang^beta+Epsilon);
                if isnan(W(i,j))
                    W(i,j)=0;
                end                
            end
        end
    else
        J=cell(N,1);
        for ii=1:N
            Xhat = X(:,neighborhood(:,ii))-repmat(X(:,ii),1,K); % shift ith pt to origin
            C = Xhat'*Xhat;                                        % local covariance
            S = diag(C);
            S = diag( 1./(1+sigmaN.*S)  ) ;
            [J{ii},~] = eigs(Xhat*(S*S')*Xhat',d);
        end
        W=zeros(N);
        for i=1:N
            ind=find(distance(i,:)>0);
            for j=ind
                [~,Angle]=prinAngles(J{i},J{j},'fast');
                Ang=sum(abs(Angle).^2);
                W(i,j)=(distance(i,j)+Epsilon)/(Ang^beta+Epsilon);
                if isnan(W(i,j))
                    W(i,j)=0;
                end                
            end
        end
    end
end

% W(W==0)=inf;