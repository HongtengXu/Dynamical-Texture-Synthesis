function out=ManiDynTexSyn(in1,inT,pr,wr,step,T,K,d,Epsilon,sigmaN,method,beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Manifold Guided Dynamic Texture Synthesis
% Inputs:
%   in1, inT: two source images
%   pr: the radius of patch
%   step: the sampling step
%   T: the number of frame to synthesize
%   K: the number of neighbors find in for each temporal patch
%   lambda: the weight of gradient domain in patch matching
%
% Output:
%   out: the output sequence
%
% The details of algorithm
% 1) Apply optical flow for patch matching
% 2) To the patches whose amplitude of flow are smaller than threshold,
% apply linear synthesize
% 3) To the patches with large optical flow, estimate the matching regions
% firstly.
% 4) Sampling regions to get k neighbors for each temporal patch.
% 5) Calculate to each neighbor, calculate their graph matrix and principal
% angles
% 6) Applying Dijkstra's minimum distance algorithm for graphs.
% 7) Stitching patches by graph-cut.
%
% Hongteng Xu
% School of ECE, Georgia Tech
% 09/20/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[R,C,L]=size(in1);
R=floor((R-2*pr-1)/step)*step+2*pr+1;
C=floor((C-2*pr-1)/step)*step+2*pr+1;
in1=in1(1:R,1:C,:);
inT=inT(1:R,1:C,:);
out=zeros(R,C,L,T);
out(:,:,:,1)=in1;
out(:,:,:,T)=inT;
[cols,rows]=meshgrid(1+pr:step:C-pr,1+pr:step:R-pr);
Mask=zeros(size(cols));
coord=[rows(:),cols(:)];
SE = strel('square', 3);

res=sum(abs(in1-inT),3);
err=[];
for i=1:length(rows(:))
    err=[err,res(coord(i,1),coord(i,2))];
end
[~,index]=sort(err,'ascend');
r=coord(index(1),1);
c=coord(index(1),2);
Mask=double(cols==c & rows==r);

patch1=in1(r-pr:r+pr,c-pr:c+pr,:);
patchT=inT(r-pr:r+pr,c-pr:c+pr,:);
for t=2:T-1               
    out(r-pr:r+pr,c-pr:c+pr,:,t)= ( (T-1-(t-1))*patch1+(t-1)*patchT )./(T-1); 
end   
MaskFormer=Mask;
NUM=1;

tic;
while sum(MaskFormer(:))<length(rows(:))
    
    Mask=imdilate(Mask,SE);
    tmp=Mask-MaskFormer;
    index=find(tmp(:)==1);
    for j=1:length(index)
        i=index(j);
        r=coord(i,1);
        c=coord(i,2);
        
        % local region estimation and create candidates
        patch1=in1(r-pr:r+pr,c-pr:c+pr,:);
        patchT=inT(r-pr:r+pr,c-pr:c+pr,:);
        
        OverLap=double(sum(out(r-pr:r+pr,c-pr:c+pr,:,2),3)>0);
        
        row1=max([r-wr,1]);
        row2=min([row1+2*wr,R]);
        col1=max([c-wr,1]);
        col2=min([col1+2*wr,C]);
        img1=in1(row1:row2,col1:col2,:);
        imgT=inT(row1:row2,col1:col2,:);
        [S1,samples1]=SampleEstimate(patch1, img1);
        [ST,samplesT]=SampleEstimate(patchT, imgT);
        Samples=[samples1,samplesT];
    
            
        % candidates selection
%       tic;
        overlap=[];

        for t=2:T-1
            tmp=im2double(out(r-pr:r+pr,c-pr:c+pr,:,t));
            overlap=[overlap,tmp(:)];
            mask = double(tmp(:)~=0);
        end
        Candidates=CandidateSelection(S1,ST,Samples,overlap,mask,K,T,0);


        [W,X,MASK]=GraphCreatation(Candidates, sigmaN, K, d, Epsilon,method,beta);


        path=PathSelection(W);


        for t=2:T-1
            patch=reshape( X(:,path(t)), [2*pr+1,2*pr+1,L] );
            pover=reshape( overlap(:,t-1), [2*pr+1,2*pr+1,L] );                 
            out(r-pr:r+pr,c-pr:c+pr,:,t)=PatchFitting(pover, patch, step, pr, L, OverLap);
        end
        NUM=NUM+1;
%                  imshow(out(:,:,:,2))
        fprintf('%d/%d,time=%d\n',NUM,size(coord,1),toc);
    end
        
    
    MaskFormer=Mask;
end
