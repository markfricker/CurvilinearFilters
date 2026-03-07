function [imOut] = BowlerHatFilter(imIn,minScale,nScales,nOrientations)
mn = min(imIn,[],'all');
mx = max(imIn,[],'all');
r = minScale:minScale+nScales; %radius of the disk
l = r.*2+1; % length of the line
o = 0:180/nOrientations:180-180/nOrientations; %number of orientation
imIn = im2single(imIn);
imol = zeros(size(imIn,1),size(imIn,2),length(r),length(o),'single');
imod = zeros(size(imIn,1),size(imIn,2),length(r),'single');
for iS=1:length(r)
    for iO=1:length(o)
        se = strel('line',l(iS),o(iO));
        imol(:,:,iS,iO) = imopen(imIn,se);
    end
    se = strel('disk',r(iS));% original
    imod(:,:,iS) = imopen(imIn,se);
end
%% Diff
imd = zeros(size(imIn,1),size(imIn,2),length(r),'single');
imm = zeros(size(imIn,1),size(imIn,2),length(r),'single');
for iS=1:length(r)
    imm(:,:,iS) = max(squeeze(imol(:,:,iS,:)),[],3);   % Max for all lines
    imd(:,:,iS) = imm(:,:,iS) - imod(:,:,iS);           % Diff betwen disk and line
end
imda = max(imd,[],3);
imOut = cast(rescale(imda,mn,mx),'like',imIn);
end
