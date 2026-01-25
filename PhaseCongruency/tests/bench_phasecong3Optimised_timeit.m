function bench_phasecong3Optimised_timeit()

    if exist('reticulate_network_test.png','file')
        im = imread('reticulate_network_test.png');
    else
        error('reticulate_network_test.png not found in current folder.');
    end
    if ndims(im)==3, im = rgb2gray(im); end

    imD = double(im);
    imS = single(imD);

    % Warm-up
    phasecong3Optimised(imD, 'precision','double', 'storeEO',true, 'storePC',false);

    cases = {
        'double, EO=true, PC=false (recommended)', @() phasecong3Optimised(imD, 'precision','double', 'storeEO',true,  'storePC',false)
        'single, EO=true, PC=false (recommended)', @() phasecong3Optimised(imS, 'precision','single', 'storeEO',true,  'storePC',false)

        'double, EO=false, PC=false',              @() phasecong3Optimised(imD, 'precision','double', 'storeEO',false, 'storePC',false)
        'single, EO=false, PC=false',              @() phasecong3Optimised(imS, 'precision','single', 'storeEO',false, 'storePC',false)

        'double, EO=true, PC=true',                @() phasecong3Optimised(imD, 'precision','double', 'storeEO',true,  'storePC',true)
        'single, EO=true, PC=true',                @() phasecong3Optimised(imS, 'precision','single', 'storeEO',true,  'storePC',true)

        % Including noiseMask output (slightly more work)
        'double + noiseMask',                      @() call_with_noisemask(imD,'double')
        'single + noiseMask',                      @() call_with_noisemask(imS,'single')
    };

    t = zeros(size(cases,1),1);

    fprintf('\nRunning timeit() benchmarks for phasecong3Optimised...\n');
    for i = 1:size(cases,1)
        t(i) = timeit(cases{i,2});
        fprintf('%2d) %-45s  %.4f s\n', i, cases{i,1}, t(i));
    end

    T = table(string(cases(:,1)), t, 'VariableNames', {'Case','Seconds'});
    disp(T);

    fprintf('\nTo clear the persistent cache: clear phasecong3Optimised\n');
end

function call_with_noisemask(im, prec)
    [~,~,~,~,~,~,~,~,~] = phasecong3Optimised(im, 'precision',prec, 'storeEO',true, 'storePC',false);
end
