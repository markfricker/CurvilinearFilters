function Jnodes = detectJunctions(response, direction, varargin)
%DETECTJUNCTIONS Detect junction nodes in curvilinear networks
%
% Returns ONE pixel per junction (node-level detector)

p = inputParser;
addParameter(p,'Radius',3);
addParameter(p,'ThresholdFrac',0.2);
addParameter(p,'AngleBins',6);
addParameter(p,'MinVotes',2);
parse(p,varargin{:});

r   = p.Results.Radius;
thr = p.Results.ThresholdFrac * max(response(:));
nb  = p.Results.AngleBins;
mv  = p.Results.MinVotes;

% ---- Step 1: pixel-level candidates (orientation voting) ----
Jpix = false(size(response));
angEdges = linspace(0,pi,nb+1);

for i = 1+r:size(response,1)-r
    for j = 1+r:size(response,2)-r
        winR = response(i-r:i+r, j-r:j+r);
        mask = winR > thr;

        if nnz(mask) < 5
            continue;
        end

        winD = direction(i-r:i+r, j-r:j+r);
        ang  = mod(winD(mask), pi);
        counts = histcounts(ang, angEdges);

        if nnz(counts > 0) >= mv
            Jpix(i,j) = true;
        end
    end
end

% ---- Step 2: cluster pixels ----
CC = bwconncomp(Jpix);

% ---- Step 3: reduce each cluster to ONE node ----
Jnodes = false(size(response));

stats = regionprops(CC,'Centroid');

for k = 1:CC.NumObjects
    c = round(stats(k).Centroid);
    % Guard bounds
    c(1) = min(max(c(1),1),size(Jnodes,2));
    c(2) = min(max(c(2),1),size(Jnodes,1));
    Jnodes(c(2),c(1)) = true;
end
end




