function tests = test_neuriteness_continuity
tests = functiontests(localfunctions);
end

function testContinuity(testCase)
I = generateTestLineImage(128,30,2);

[N,~] = neuriteness2D(I,2);

bw = N > 0.4*max(N(:));
cc = bwconncomp(bw);

verifyLessThanOrEqual(testCase,cc.NumObjects,2);
end
