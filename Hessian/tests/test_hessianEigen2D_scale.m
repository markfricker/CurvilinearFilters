function tests = test_hessianEigen2D_scale
tests = functiontests(localfunctions);
end

function testScaleNormalization(testCase)
I = generateTestBlobImage(64,4);

[L1a,~,~,~] = hessianEigen2D(I,2);
[L1b,~,~,~] = hessianEigen2D(I,4);

ratio = max(abs(L1a(:))) / max(abs(L1b(:)));

verifyGreaterThan(testCase, ratio, 0.2);
verifyLessThan(testCase, ratio, 5);
end
