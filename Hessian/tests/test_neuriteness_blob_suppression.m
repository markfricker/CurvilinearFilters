function tests = test_neuriteness_blob_suppression
tests = functiontests(localfunctions);
end

function testBlobSuppression(testCase)
Iline = generateTestLineImage(128,0,2);
Iblob = generateTestBlobImage(128,6);

Nline = neuriteness2D(Iline,2);
Nblob = neuriteness2D(Iblob,2);

lineSupport = nnz(Nline > 0.5);
blobSupport = nnz(Nblob > 0.5);

verifyGreaterThan(testCase, lineSupport, 3 * blobSupport);
end
