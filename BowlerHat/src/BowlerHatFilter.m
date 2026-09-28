function [imOut] = BowlerHatFilter(imIn,minScale,nScales,nOrientations)
% BowlerHatFilter  Multiscale bowler-hat transform for curvilinear enhancement.
%
%   imOut = BowlerHatFilter(imIn, minScale, nScales, nOrientations)
%   imOut = BowlerHatFilter(imIn, radii, nOrientations)
%
% OVERVIEW
%   The bowler-hat transform enhances thin curvilinear structures (vessels,
%   tubules, fibres) using two banks of morphological openings: a bank of
%   disk structuring elements (which suppress structures narrower than the
%   disk regardless of orientation) and a bank of oriented line structuring
%   elements (which preserve structures aligned with the line but suppress
%   those that are not). At each scale, the response is the difference
%   between the maximum line-opening (over all orientations) and the
%   disk-opening; taking the maximum of this difference across scales
%   yields a near-isotropic ridge/line map that is robust to junctions and
%   branching, where Hessian-based filters (e.g. Frangi vesselness) tend to
%   under-respond.
%
% INPUTS
%   imIn          - 2D grayscale image.
%   minScale      - Smallest disk/line radius (px).
%   nScales       - Number of additional scales beyond minScale (radii run
%                   minScale:minScale+nScales).
%   nOrientations - Number of line orientations sampled over [0,180) deg.
%   radii         - (3-argument form) explicit vector of disk/line radii
%                   (px), e.g. sigmaMin:sigmaStep:sigmaMax. Rounded to
%                   whole pixels (strel('disk') needs an integer radius)
%                   and de-duplicated; values below 1 are raised to 1.
%
% OUTPUT
%   imOut - Bowler-hat response, rescaled back to the input's intensity
%           range and cast to the input's class.
%
% REFERENCES
%   Sazak, C., Nelson, C. J., & Obara, B. (2019). The multiscale bowler-hat
%   transform for blood vessel enhancement in retinal images. Pattern
%   Recognition, 88, 739-750. https://doi.org/10.1016/j.patcog.2018.10.011
%
mn = min(imIn,[],'all');
mx = max(imIn,[],'all');
if nargin == 3
    % BowlerHatFilter(imIn, radii, nOrientations)
    nOrientations = nScales;
    r = unique(max(1, round(double(minScale(:)'))));
else
    r = minScale:minScale+nScales; %radius of the disk
end
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
