classdef TestBowlerHat < matlab.unittest.TestCase
% TestBowlerHat  Unit tests for BowlerHatFilter.
%
% BowlerHatFilter detects elongated curvilinear structures (rods/tubules)
% by computing the multi-scale difference between the maximum line-opening
% and the disk-opening of an image.
%
% USAGE
%   results = runtests('BowlerHat/tests/TestBowlerHat');
%   table(results)
%
% REQUIREMENTS
%   Image Processing Toolbox (imopen, strel, imgaussfilt)

    % =====================================================================
    % Path setup
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            bhRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(bhRoot, 'src'));
        end
    end

    % =====================================================================
    % Shared test images
    % =====================================================================
    methods (Static, Access = private)

        function I = rodImage()
            % Short bright horizontal rod (16 px) on dark background.
            % Gaussian cross-section (sigma 3 px), centred at row 32.
            [xx, yy] = meshgrid(1:64, 1:64);
            I = single(exp(-(yy - 32).^2 / (2*3^2)) .* (xx >= 24 & xx <= 39));
            I = I / max(I(:));
        end

        function I = blobImage()
            % Isotropic Gaussian blob – no preferred orientation.
            [xx, yy] = meshgrid(1:64, 1:64);
            I = single(exp(-((xx-32).^2 + (yy-32).^2) / (2*6^2)));
            I = I / max(I(:));
        end

    end

    % =====================================================================
    % Default test parameters (small bank for speed)
    % =====================================================================
    properties (Constant)
        MinScale    = 2
        NScales     = 2
        NOrient     = 4
    end

    % =====================================================================
    % Tests
    % =====================================================================
    methods (Test)

        function testBH_smoke(tc)
            % No-crash smoke test
            out = BowlerHatFilter( ...
                TestBowlerHat.rodImage(), ...
                tc.MinScale, tc.NScales, tc.NOrient);
            tc.verifyNotEmpty(out);
        end

        function testBH_outputSize(tc)
            I   = TestBowlerHat.rodImage();
            out = BowlerHatFilter(I, tc.MinScale, tc.NScales, tc.NOrient);
            tc.verifySize(out, size(I));
        end

        function testBH_outputClass_single(tc)
            % Single input → single output (cast(...,'like',imIn))
            I   = TestBowlerHat.rodImage();
            out = BowlerHatFilter(I, tc.MinScale, tc.NScales, tc.NOrient);
            tc.verifyClass(out, 'single');
        end

        function testBH_outputClass_uint8(tc)
            % uint8 input → uint8 output
            I   = uint8(255 * TestBowlerHat.rodImage());
            out = BowlerHatFilter(I, tc.MinScale, tc.NScales, tc.NOrient);
            tc.verifyClass(out, 'uint8');
        end

        function testBH_outputRangeWithinInput(tc)
            % Output is rescaled to [min(input), max(input)]
            I   = TestBowlerHat.rodImage();
            out = BowlerHatFilter(I, tc.MinScale, tc.NScales, tc.NOrient);
            tc.verifyGreaterThanOrEqual(double(min(out(:))), double(min(I(:))) - 1e-4);
            tc.verifyLessThanOrEqual(   double(max(out(:))), double(max(I(:))) + 1e-4);
        end

        function testBH_rodResponse(tc)
            % The filter should respond more strongly along the rod than
            % in a flat background region
            I   = TestBowlerHat.rodImage();
            out = BowlerHatFilter(I, tc.MinScale, tc.NScales, tc.NOrient);
            rodMean = mean(out(30:34, 24:39), 'all');
            bgMean  = mean(out(1:10,  1:10 ), 'all');
            tc.verifyGreaterThan(double(rodMean), double(bgMean));
        end

        function testBH_rodVsBlob(tc)
            % BowlerHat is designed for elongated structures; it should
            % give a stronger max response on a rod than on an isotropic blob
            Irod  = TestBowlerHat.rodImage();
            Iblob = TestBowlerHat.blobImage();
            outRod  = BowlerHatFilter(Irod,  tc.MinScale, tc.NScales, tc.NOrient);
            outBlob = BowlerHatFilter(Iblob, tc.MinScale, tc.NScales, tc.NOrient);
            tc.verifyGreaterThan(double(max(outRod(:))), double(max(outBlob(:))));
        end

    end

end
