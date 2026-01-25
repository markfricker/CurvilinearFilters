function tests = test_neuriteness_direction
tests = functiontests(localfunctions);
end

function testDirection(testCase)
I = generateTestLineImage(128,0,2);

[~,D] = neuriteness2D(I,2);
mask = abs(D) < pi/4;

verifyGreaterThan(testCase,nnz(mask),500);
end
