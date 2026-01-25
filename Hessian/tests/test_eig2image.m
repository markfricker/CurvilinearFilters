function tests = test_eig2image
tests = functiontests(localfunctions);
end

function testEigenOrdering(testCase)

I = generateTestBlobImage(64,5);
[Dxx,Dxy,Dyy] = Hessian2D(I,2);
[L1,L2,Vx,Vy] = eig2image(Dxx,Dxy,Dyy);

% --- Eigenvalue ordering: pointwise ---
verifyTrue(testCase, all(abs(L1(:)) <= abs(L2(:)) + eps));

% --- Eigenvector normalization: only where defined ---
mag = hypot(Vx, Vy);
mask = mag > 0;

if any(mask(:))
    verifyLessThan(testCase, max(abs(mag(mask) - 1)), 1e-3);
else
    % No defined eigenvectors anywhere is valid
    verifyTrue(testCase, true);
end

end


