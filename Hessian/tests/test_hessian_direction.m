function tests = test_hessian_direction
tests = functiontests(localfunctions);
end

function testDirectionAlongLine(testCase)
I = generateTestLineImage(128,0,2);

[V,~,D] = hessian2DFilters(I, ...
    'FilterType','vesselness', ...
    'Sigmas',2, ...
    'Parameters',struct('beta',0.5,'c',15));

mask = V > 0.5*max(V(:));
meanAngle = atan2(mean(sin(D(mask))),mean(cos(D(mask))));

verifyLessThan(testCase,abs(meanAngle),pi/8);
end
