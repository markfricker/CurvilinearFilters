function tests = test_Hessian2D
tests = functiontests(localfunctions);
end

function testConstantImage(testCase)
I = ones(64);
[Dxx,Dxy,Dyy] = Hessian2D(I,2);

energy = mean(abs([Dxx(:);Dxy(:);Dyy(:)]));
verifyLessThan(testCase, energy, 0.1);
end
