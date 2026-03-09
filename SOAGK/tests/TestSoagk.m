classdef TestSoagk < matlab.unittest.TestCase
% TestSoagk  Unit tests for AGlineDetectorSteerable and AGlineDetectorSteerableConv2.
%
% Migrated from functiontests pattern to matlab.unittest.TestCase class.
%
% USAGE
%   results = runtests('SOAGK/tests/TestSoagk');
%   table(results)

    % =====================================================================
    % Path setup
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            soagkRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(soagkRoot, 'src'));
            addpath(fullfile(soagkRoot, 'utils'));
        end
    end

    % =====================================================================
    % Test image helper (static private)
    % =====================================================================
    methods (Static, Access = private)
        function I = makeHLine(sz, rowPos, width)
            % Normalised image with a horizontal bar of given width.
            I  = zeros(sz, sz, 'single');
            r0 = max(1,  rowPos - floor(width/2));
            r1 = min(sz, rowPos + floor(width/2));
            I(r0:r1, :) = 1;
            I = imgaussfilt(I, 0.8);
            I = I ./ max(I(:));
        end
    end

    % =====================================================================
    % Tests
    % =====================================================================
    methods (Test)

        function testOutputSize(tc)
            I      = TestSoagk.makeHLine(64, 32, 2);
            sigmas = 2;
            thetas = 0 : pi/4 : pi - pi/4;
            rhos   = 2;
            [lineMap, dirMap] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos);
            tc.verifySize(lineMap, size(I));
            tc.verifySize(dirMap,  size(I));
        end

        function testLineMapRange(tc)
            I      = TestSoagk.makeHLine(64, 32, 2);
            sigmas = 2;
            thetas = 0 : pi/4 : pi - pi/4;
            rhos   = 2;
            [lineMap, ~] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos);
            tc.verifyGreaterThanOrEqual(min(lineMap(:)), 0);
            tc.verifyLessThanOrEqual(   max(lineMap(:)), 1 + 1e-6);
        end

        function testHorizontalLineResponse(tc)
            I      = TestSoagk.makeHLine(128, 64, 3);
            sigmas = 1:3;
            thetas = 0 : pi/6 : pi - pi/6;
            rhos   = [1 2];
            [lineMap, ~] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos);
            lineResp = mean(lineMap(60:68, 30:100), 'all');
            bgResp   = mean(lineMap(1:20,  1:30),   'all');
            tc.verifyGreaterThan(lineResp, 5 * bgResp + 1e-6);
        end

        function testOrientationAccuracy(tc)
            % Horizontal line → dirMap at line pixels should be near pi/2
            % (SOAGK maximises the normal-direction second derivative)
            I      = TestSoagk.makeHLine(128, 64, 3);
            sigmas = 2;
            thetas = 0 : pi/12 : pi - pi/12;
            rhos   = 2;
            [lineMap, dirMap] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos);
            mask        = lineMap > 0.5;
            residual    = abs(dirMap(mask) - pi/2);
            medResidual = median(residual);
            tc.verifyLessThan(medResidual, pi/6);   % within 30 deg of pi/2
        end

        function testConv2MatchesImfilter(tc)
            % Interior of conv2 and imfilter variants should agree
            I      = TestSoagk.makeHLine(64, 32, 2);
            sigmas = 1:2;
            thetas = 0 : pi/4 : pi - pi/4;
            rhos   = [1 2];
            [L1, D1] = AGlineDetectorSteerable(I,        sigmas, thetas, rhos, 'single');
            [L2, D2] = AGlineDetectorSteerableConv2(I,   sigmas, thetas, rhos, 'single');
            border = 12;
            r = border+1 : size(I,1)-border;
            c = border+1 : size(I,2)-border;
            tc.verifyEqual(L1(r,c), L2(r,c), 'AbsTol', single(1e-3));
            tc.verifyEqual(D1(r,c), D2(r,c), 'AbsTol', single(1e-3));
        end

        function testSingleVsDoubleClose(tc)
            I      = TestSoagk.makeHLine(64, 32, 2);
            sigmas = 2;
            thetas = 0 : pi/4 : pi - pi/4;
            rhos   = 2;
            [Ls, ~] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos, 'single');
            [Ld, ~] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos, 'double');
            tc.verifyEqual(double(Ls), Ld, 'AbsTol', 1e-3);
        end

    end

end
