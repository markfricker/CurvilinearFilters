function tests = test_vesselness_polarity
tests = functiontests(localfunctions);
end

function testPolarity(testCase)
I = generateTestLineImage(128,0,2);

[V1,~,~] = hessian2DFilters(I,'FilterType','vesselness','Sigmas',2,'WhiteOnDark',true);
[V2,~,~] = hessian2DFilters(I,'FilterType','vesselness','Sigmas',2,'WhiteOnDark',false);

verifyGreaterThan(testCase,max(V1(:)),5*max(V2(:)));
end
d