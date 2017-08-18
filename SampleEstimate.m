function [target,samples]=SampleEstimate(patch, img)



ps=size(patch,1);
samples=[];
for l=1:size(img,3)
    samples=[samples;im2col(img(:,:,l),[ps,ps],'sliding')];
end
target=patch(:);

res=sum(abs(samples-repmat(target,[1,size(samples,2)])));
ind= res~=0;
samples=samples(:,ind);

