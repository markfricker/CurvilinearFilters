function tests = test_vesselness_basic
tests = functiontests(localfunctions);
end

function testLineDetection(testCase)
I = generateTestLineImage(128,0,2);

[V,~,~] = hessian2DFilters(I, ...
    'FilterType','vesselness', ...
    'Sigmas',1:4, ...
    'Parameters',struct('beta',0.5,'c',15));

ridgeVal = mean(V(50:80,50:80),'all');
bgVal    = mean(V(1:20,1:20),'all');

verifyGreaterThan(testCase,ridgeVal,5*bgVal);
end
