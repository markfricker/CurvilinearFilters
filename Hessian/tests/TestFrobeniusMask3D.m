classdef TestFrobeniusMask3D < matlab.unittest.TestCase
% TestFrobeniusMask3D  Unit tests for localFrobeniusMask3D.
%
% USAGE
%   results = runtests('Hessian/tests/TestFrobeniusMask3D');
%
% REQUIREMENTS
%   Image Processing Toolbox (imgaussfilt3)
%   Segmentation_sandbox/src/ on the path (globalThresholdFast)

    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            hRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(hRoot, 'src'));
            addpath(fullfile(hRoot, 'src', 'engine'));
            segRoot = fullfile(hRoot, '..', '..', 'Segmentation_sandbox', 'src');
            tc.assertTrue(isfolder(segRoot), sprintf('Segmentation_sandbox src not found at %s', segRoot));
            addpath(segRoot);
            tc.assertTrue(exist('globalThresholdFast', 'file') > 0, ...
                'globalThresholdFast not resolvable after path setup.');
        end
    end

    methods (Test)

        function testGate_zerosBackground_keepsStructure(tc)
            [im, tubeMask] = tc.tubeWithBackground(48, 3);
            sigmas = 1:4;

            imf = ones(size(im), 'single');   % isolate the gate's own effect
            imfMasked = localFrobeniusMask3D(imf, im, sigmas, 2);

            % Deep background, far from the tube, must be zeroed.
            bgRegion = false(size(im));
            bgRegion(1:6, 1:6, 1:6) = true;   % a corner, nowhere near the centred tube
            tc.verifyEqual(nnz(imfMasked(bgRegion)), 0, ...
                'far-background voxels must be hard-zeroed by the gate');

            % Tube core must survive.
            [nY, nX, nZ] = size(im);
            core = false(size(im));
            core(round(nY/2), round(nX/2), round(nZ/2)) = true;
            tc.verifyGreaterThan(imfMasked(core), 0, ...
                'the tube core must survive the gate');
        end

        function testFrobDivision_isMonotonicallyPermissive(tc)
            [im, ~] = tc.tubeWithBackground(48, 3);
            sigmas = 1:4;
            imf = ones(size(im), 'single');

            nSurvive = zeros(1,3);
            divisions = [1, 2, 4];
            for k = 1:numel(divisions)
                masked = localFrobeniusMask3D(imf, im, sigmas, divisions(k));
                nSurvive(k) = nnz(masked);
            end

            tc.verifyTrue(all(diff(nSurvive) >= 0), ...
                'larger frobDivision (more permissive) must never let FEWER voxels survive');
            tc.verifyGreaterThan(nSurvive(end), nSurvive(1), ...
                'frobDivision=4 must let strictly more voxels through than frobDivision=1 on this phantom');
        end

        function testDefaultFrobDivision_isTwo(tc)
            [im, ~] = tc.tubeWithBackground(32, 3);
            sigmas = 1:3;
            imf = ones(size(im), 'single');

            explicitDefault = localFrobeniusMask3D(imf, im, sigmas, 2);
            implicitDefault = localFrobeniusMask3D(imf, im, sigmas);

            tc.verifyEqual(implicitDefault, explicitDefault, ...
                'omitting frobDivision must behave exactly like frobDivision=2');
        end

        function testOutputPreservesRealResponseValues(tc)
            % Not just a 0/1 mask applied to a constant -- confirm real
            % (non-uniform) enhance values pass through unchanged where kept.
            [im, ~] = tc.tubeWithBackground(32, 3);
            sigmas = 1:3;
            rng(1);
            imf = single(rand(size(im)));

            masked = localFrobeniusMask3D(imf, im, sigmas, 2);
            keptMask = (masked ~= 0);
            tc.verifyEqual(masked(keptMask), imf(keptMask), ...
                'surviving voxels must retain their original imf value, unmodified');
        end
    end

    methods (Access = private, Static)
        function [V, tubeMask] = tubeWithBackground(sz, radius)
            [X, Y, ~] = ndgrid(1:sz, 1:sz, 1:sz);
            cx = (sz+1)/2; cy = (sz+1)/2;
            tubeMask = hypot(X-cx, Y-cy) <= radius;
            V = single(tubeMask);
            V = imgaussfilt3(V, 1);
        end
    end
end
