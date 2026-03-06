function tests = test_AGlineDetector
% test_AGlineDetector
% Unit tests for AGlineDetectorSteerable and AGlineDetectorSteerableConv2.
%
% Run with: results = runtests('test_AGlineDetector')
tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------------------
% Setup: add src to path
% ---------------------------------------------------------------------------
function setupOnce(testCase) %#ok<INUSD>
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir,'..','src'));
addpath(fullfile(thisDir,'..','utils'));
end

% ---------------------------------------------------------------------------
% Tests
% ---------------------------------------------------------------------------

function testOutputSize(testCase)
% lineMap and dirMap must be the same size as the input image.
I      = makeHLine(64, 32, 2);
sigmas = 2;
thetas = 0 : pi/4 : pi - pi/4;
rhos   = 2;

[lineMap, dirMap] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos);

verifyEqual(testCase, size(lineMap), size(I));
verifyEqual(testCase, size(dirMap),  size(I));
end

function testLineMapRange(testCase)
% lineMap output must lie in [0, 1].
I      = makeHLine(64, 32, 2);
sigmas = 2;
thetas = 0 : pi/4 : pi - pi/4;
rhos   = 2;

[lineMap, ~] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos);

verifyGreaterThanOrEqual(testCase, min(lineMap(:)), 0);
verifyLessThanOrEqual(testCase,    max(lineMap(:)), 1 + 1e-6);
end

function testHorizontalLineResponse(testCase)
% A horizontal line should produce a strong response at the line pixels
% and a weak response in the background.
I = makeHLine(128, 64, 3);   % horizontal bar at row 64

sigmas = 1:3;
thetas = 0 : pi/6 : pi - pi/6;
rhos   = [1 2];

[lineMap, ~] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos);

lineResp = mean(lineMap(60:68, 30:100), 'all');
bgResp   = mean(lineMap(1:20,  1:30),   'all');

verifyGreaterThan(testCase, lineResp, 5 * bgResp + 1e-6);
end

function testOrientationAccuracy(testCase)
% For a horizontal line, the SOAGK detector maximises the negative second
% derivative in the direction NORMAL to the line (theta = pi/2), so
% dirMap should be close to pi/2 at line pixels.
I      = makeHLine(128, 64, 3);
sigmas = 2;
thetas = 0 : pi/12 : pi - pi/12;
rhos   = 2;

[lineMap, dirMap] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos);

% Identify strong-response pixels along the line
mask = lineMap > 0.5;
linePixelDirs = dirMap(mask);

% dirMap stores the theta of max second-derivative response (normal to line)
residual    = abs(linePixelDirs - pi/2);
medResidual = median(residual);

verifyLessThan(testCase, medResidual, pi/6);   % within 30 degrees of pi/2
end

function testConv2MatchesImfilter(testCase)
% The conv2 and imfilter implementations agree in the image interior.
% They differ at borders: imfilter uses replicate padding, conv2 uses
% zero-padding via 'same' mode.  Exclude a conservative border.
I      = makeHLine(64, 32, 2);
sigmas = 1:2;
thetas = 0 : pi/4 : pi - pi/4;
rhos   = [1 2];

[L1, D1] = AGlineDetectorSteerable(I, sigmas, thetas, rhos, 'single');
[L2, D2] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos, 'single');

% Max kernel half-width for sigma=2, rho=1: floor(round(9*2)/2) = 9
border = 12;
r = border+1 : size(I,1)-border;
c = border+1 : size(I,2)-border;

verifyEqual(testCase, L1(r,c), L2(r,c), 'AbsTol', single(1e-3), ...
    'conv2 and imfilter lineMap interior should agree');
verifyEqual(testCase, D1(r,c), D2(r,c), 'AbsTol', single(1e-3), ...
    'conv2 and imfilter dirMap interior should agree');
end

function testSingleVsDoubleClose(testCase)
% Single and double precision should give numerically close results.
I      = makeHLine(64, 32, 2);
sigmas = 2;
thetas = 0 : pi/4 : pi - pi/4;
rhos   = 2;

[Ls, ~] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos, 'single');
[Ld, ~] = AGlineDetectorSteerableConv2(I, sigmas, thetas, rhos, 'double');

verifyEqual(testCase, double(Ls), Ld, 'AbsTol', 1e-3);
end

% ---------------------------------------------------------------------------
% Local helpers
% ---------------------------------------------------------------------------

function I = makeHLine(sz, rowPos, width)
% Create a normalized image with a horizontal bar of given width.
I = zeros(sz, sz, 'single');
r0 = rowPos - floor(width/2);
r1 = rowPos + floor(width/2);
r0 = max(1, r0);
r1 = min(sz, r1);
I(r0:r1, :) = 1;
I = imgaussfilt(I, 0.8);   % slight blur to avoid Gibbs ringing
I = I ./ max(I(:));
end
