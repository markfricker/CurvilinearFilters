classdef TestSteerGauss < matlab.unittest.TestCase
% TestSteerGauss  Unit tests for steerGaussEnhance.
%
% The 'sampled' mode is checked against a self-contained brute-force
% reference that reproduces the original per-orientation implementation
% (three full 2-D kernel correlations per (scale, theta), as in
% steerGaussFilterOrder2 from MycNetAnalysis) -- so the regression check
% does not need that third-party project on the path.
%
% USAGE
%   results = runtests('SOAGK/tests/TestSteerGauss');
%   table(results)

    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            soagkRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(soagkRoot, 'src'));
        end
    end

    methods (Static, Access = private)
        function I = makeTestImage(sz)
            rng(3);
            I = zeros(sz, 'single');
            I(round(sz(1)/2)-1:round(sz(1)/2)+1, :) = 1;        % horizontal bar
            I(:, round(sz(2)/3)) = 1;                           % vertical line
            I = imgaussfilt(I, 1) + 0.05*randn(sz, 'single');
            I = mat2gray(I);
        end

        function [L, D] = bruteForce(I, sigmas, thetas)
            % Original algorithm: 2-D kernels per (scale, theta), stacked max.
            I = double(I);
            nS = numel(sigmas); nT = numel(thetas);
            respAll = zeros([size(I), nS*nT], 'single');
            idx = 0;
            for is = 1:nS
                s  = sigmas(is);
                Wx = max(1, floor(4*s));
                [xx, yy] = meshgrid(-Wx:Wx, -Wx:Wx);
                g0  = exp(-(xx.^2+yy.^2)/(2*s^2))/(s*sqrt(2*pi));
                G2a = -g0/s^2 + g0.*xx.^2/s^4;
                G2b =  g0.*xx.*yy/s^4;
                G2c = -g0/s^2 + g0.*yy.^2/s^4;
                I2a = imfilter(I, G2a, 'same', 'replicate');
                I2b = imfilter(I, G2b, 'same', 'replicate');
                I2c = imfilter(I, G2c, 'same', 'replicate');
                for it = 1:nT
                    idx = idx + 1;
                    t = -thetas(it)*pi/180;
                    J = cos(t)^2*I2a + sin(t)^2*I2c - 2*cos(t)*sin(t)*I2b;
                    respAll(:,:,idx) = -single(J);
                end
            end
            [L, maxIdx] = max(respAll, [], 3);
            D = single(thetas(mod(maxIdx - 1, nT) + 1));
            L(L < 0) = 0;
            L = L ./ max(L(:));
        end
    end

    methods (Test)

        function sampledMatchesOriginalAlgorithm(tc)
            I      = TestSteerGauss.makeTestImage([96 128]);
            sigmas = 0.5:0.5:3;
            thetas = 0:15:345;
            [Lref, Dref] = TestSteerGauss.bruteForce(I, sigmas, thetas);
            [L, D] = steerGaussEnhance(I, sigmas, thetas, 'single');
            tc.verifyClass(L, 'single');
            tc.verifyEqual(L, Lref, 'AbsTol', single(1e-5));
            tc.verifyGreaterThanOrEqual(mean(D(:) == Dref(:)), 0.999);
        end

        function defaultModeIsSampled(tc)
            I = TestSteerGauss.makeTestImage([64 64]);
            [L1, D1] = steerGaussEnhance(I, [1 2], 0:15:345, 'single');
            [L2, D2] = steerGaussEnhance(I, [1 2], 0:15:345, 'single', 'sampled');
            tc.verifyEqual(L1, L2);
            tc.verifyEqual(D1, D2);
        end

        function analyticMatchesDenseSampling(tc)
            I      = TestSteerGauss.makeTestImage([96 128]);
            sigmas = [1 2 3];
            [La, Da] = steerGaussEnhance(I, sigmas, [], 'single', 'analytic');
            [Ld, Dd] = steerGaussEnhance(I, sigmas, 0:0.5:179.5, 'single');
            tc.verifyEqual(La, Ld, 'AbsTol', single(1e-3));
            strong = Ld > 0.2;
            dA = abs(mod(Da(strong) - Dd(strong) + 90, 180) - 90);
            tc.verifyLessThan(max(dA), 1);
            tc.verifyGreaterThanOrEqual(min(Da(:)), single(0));
            tc.verifyLessThan(max(Da(:)), single(180));
        end

        function brightRidgePolarityAndDirection(tc)
            % Bright horizontal bar: strong response on the bar, orientation
            % of max curvature is ACROSS it (90 deg).
            I = zeros(64, 'single');
            I(31:33, :) = 1;
            I = imgaussfilt(I, 1);
            for mode = {'sampled', 'analytic'}
                [L, D] = steerGaussEnhance(I, 1.5, 0:15:345, 'single', mode{1});
                tc.verifyGreaterThan(L(32, 32), single(0.9), mode{1});
                tc.verifyLessThan(L(10, 32), single(0.05), mode{1});
                tc.verifyEqual(mod(D(32, 32), 180), single(90), 'AbsTol', single(1), mode{1});
            end
        end

        function unknownModeErrors(tc)
            tc.verifyError(@() steerGaussEnhance(zeros(8,'single'), 1, 0, 'single', 'bogus'), ...
                'steerGaussEnhance:mode');
        end
    end
end
