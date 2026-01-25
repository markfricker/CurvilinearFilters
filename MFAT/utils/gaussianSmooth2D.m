function Iblur = gaussianSmooth2D(I, sigma)
% gaussianSmooth2D  Separable gaussian smoothing wrapper
% Uses imgaussfilt if available (MATLAB), otherwise falls back to conv2.

    if sigma <= 0
        Iblur = I;
        return;
    end

    try
        % imgaussfilt supports single/double in modern MATLAB
        Iblur = imgaussfilt(I, sigma);
    catch
        % fallback separable kernel
        siz = max(3, 2*ceil(3*sigma)+1);
        if isvector(I)
            % 1D fallback (not expected)
            h = fspecial('gaussian', siz, sigma);
            Iblur = conv(I, h, 'same');
        else
            H = fspecial('gaussian', siz, sigma);
            Iblur = imfilter(I, H, 'same', 'replicate');
        end
    end
end
