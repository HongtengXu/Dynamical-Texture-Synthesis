%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test the manifold based algorithm for dynamic texture synthesis via
% extremely few samples.
%
% Hongteng Xu
% School of ECE, Georgia Tech
% 09/20/2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

T=30;
% radius of patch
pr=6;
% radius of search window
wr=pr*2;
% sampling step
step=pr+1;

% parameters
K=5;
beta=1;
d=2;
Epsilon=0.01;
sigmaN=0.01;
scale=0.5;

% methodology
method=cell(3,1);
method{1}='Eu';
method{2}='Ag';
method{3}='Cu';

for ii=1
    fn1=sprintf('0 (%d).avi',ii);
    inObj = VideoReader(fn1);
    nFrames = inObj.NumberOfFrames;
    R = inObj.Height;
    C = inObj.Width;
    
    input=cell(T,1);
    for t=1:T
        input{t}=imresize( im2double( read(inObj, t) ), scale );  
    end
    in1=input{1};
    inT=input{T};
    

    out=ManiDynTexSyn(in1,inT,pr,wr,step,T,K,d,Epsilon,sigmaN,method{3},beta);   
    filename=sprintf('MDTS_%d_%s.avi',ii,method{3});
            
    outObj = VideoWriter(filename,'Uncompressed AVI');
    outObj.FrameRate = 30;
    open(outObj);
    for t = 1 : T     
        tmpOut=out(:,:,:,t);      
        tmpIn=input{t};
        tmpIn=tmpIn(1:size(tmpOut,1),1:size(tmpOut,2),:);
        tmpIn=padarray(tmpIn,[2,2]);
        tmpOut=padarray(tmpOut,[2,2]);
        writeVideo(outObj,im2uint8([tmpIn,tmpOut]));
    end
    close(outObj);
        
end