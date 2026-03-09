classdef TestNERdy < matlab.unittest.TestCase
% TestNERdy  Unit tests for nERdyEnhance.
%
% Tests that do NOT require the nERdy+ Python environment:
%   testNERdy_badInput_3D_errors    – 3-D input throws nERdyEnhance:badInput
%   testNERdy_scriptNotFound_errors – wrong path throws nERdyEnhance:notFound
%
% Tests that DO require the nERdy+ Python environment (skipped when absent):
%   testNERdy_smoke                 – no-crash run
%   testNERdy_outputSize            – output matches input size
%   testNERdy_outputClass           – output is single
%   testNERdy_binaryOutput          – output contains only 0s and 1s
%
% USAGE
%   results = runtests('nERdy/tests/TestNERdy');
%   table(results)
%
% REQUIREMENTS (for functional tests)
%   Python 3.9+ with torch, torchvision, scikit-image, pillow, numpy
%   nERdy+ source tree with pre-trained weights

    % =====================================================================
    % Path setup
    % =====================================================================
    methods (TestClassSetup)
        function addSrcPath(tc) %#ok<MANU>
            nerdyRoot = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(nerdyRoot, 'src'));
        end
    end

    % =====================================================================
    % Helper: resolve whether nERdy+ is available
    % =====================================================================
    methods (Static, Access = private)

        function tf = isNERdyAvailable()
            % Returns true if nerdy_infer.py can be found at the default path.
            srcDir    = fileparts(which('nERdyEnhance'));
            if isempty(srcDir)
                tf = false; return;
            end
            nerdyInfer = fullfile(srcDir, '..', '..', '..', ...
                                  'third party', 'nERdy integration', 'nerdy_infer.py');
            tf = isfile(nerdyInfer);
        end

        function I = grayImage()
            % Small synthetic 2-D grayscale image.
            [x, y] = meshgrid(1:64, 1:64);
            I = single(exp(-((x-32).^2 + (y-32).^2) / (2*8^2)));
        end

    end

    % =====================================================================
    % Error-condition tests (no Python / nERdy required)
    % =====================================================================
    methods (Test)

        function testNERdy_badInput_3D_errors(tc)
            % 3-D input must throw nERdyEnhance:badInput immediately,
            % before any Python subprocess is called.
            I3D = repmat(TestNERdy.grayImage(), [1 1 3]);
            tc.verifyError(@() nERdyEnhance(I3D), 'nERdyEnhance:badInput');
        end

        function testNERdy_scriptNotFound_errors(tc)
            % Wrong nerdy_infer path must throw nERdyEnhance:notFound.
            params.nerdyInfer = fullfile(tempdir(), 'nonexistent_nerdy_infer.py');
            tc.verifyError( ...
                @() nERdyEnhance(TestNERdy.grayImage(), params), ...
                'nERdyEnhance:notFound');
        end

    end

    % =====================================================================
    % Functional tests (skipped when nERdy+ is not installed)
    % =====================================================================
    methods (Test)

        function testNERdy_smoke(tc)
            tc.assumeTrue(TestNERdy.isNERdyAvailable(), ...
                'nERdy+ (nerdy_infer.py) not found — test skipped');
            R = nERdyEnhance(TestNERdy.grayImage());
            tc.verifyNotEmpty(R);
        end

        function testNERdy_outputSize(tc)
            tc.assumeTrue(TestNERdy.isNERdyAvailable(), ...
                'nERdy+ not found — test skipped');
            I = TestNERdy.grayImage();
            R = nERdyEnhance(I);
            tc.verifySize(R, size(I));
        end

        function testNERdy_outputClass(tc)
            tc.assumeTrue(TestNERdy.isNERdyAvailable(), ...
                'nERdy+ not found — test skipped');
            R = nERdyEnhance(TestNERdy.grayImage());
            tc.verifyClass(R, 'single');
        end

        function testNERdy_binaryOutput(tc)
            % nERdyEnhance always returns a binary [0,1] map
            tc.assumeTrue(TestNERdy.isNERdyAvailable(), ...
                'nERdy+ not found — test skipped');
            R = nERdyEnhance(TestNERdy.grayImage());
            uniqueVals = unique(R(:));
            tc.verifyTrue(all(uniqueVals == 0 | uniqueVals == 1), ...
                'Output must be binary (0s and 1s only)');
        end

    end

end
