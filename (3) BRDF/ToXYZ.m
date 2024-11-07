function xyz = ToXYZ(spectrum, lambdanum, lambdastart, lambdaend, white)
    cie_y_integral = 106.856895;
    [row, ~] = size(spectrum);
    xyz = zeros(row,3);

    % compute matching function
    CIE_lambda = 360:830;
    nCIESamples = length(CIE_lambda);

    fileID = fopen('CIEX.binary');
    CIE_X = fread(fileID,'float');
    fclose(fileID);
    fileID = fopen('CIEY.binary');
    CIE_Y = fread(fileID,'float');
    fclose(fileID);
    fileID = fopen('CIEZ.binary');
    CIE_Z = fread(fileID,'float');
    fclose(fileID);

    X = zeros(lambdanum, 1);
    Y = zeros(lambdanum, 1);
    Z = zeros(lambdanum, 1);
    for i = 0:lambdanum-1
        w10 = lerp(i/lambdanum, lambdastart, lambdaend);
        w11 = lerp((i+1)/lambdanum, lambdastart, lambdaend);
        X(i+1) = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, w10, w11);
        Y(i+1) = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, w10, w11);
        Z(i+1) = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, w10, w11);
    end

    for i = 1:lambdanum
        xyz(:, 1) = xyz(:, 1) + X(i) * spectrum(:, i) * white(i);
        xyz(:, 2) = xyz(:, 2) + Y(i) * spectrum(:, i) * white(i);
        xyz(:, 3) = xyz(:, 3) + Z(i) * spectrum(:, i) * white(i);
    end
    scale = (lambdaend - lambdastart)/ (cie_y_integral * lambdanum);
    xyz = xyz * scale;
end

function average = AverageSpectrumSamples(lambda, vals, n, lambdastart, lambdaend)
    sum = 0;
    if lambdastart < lambda(1)
        sum = sum + vals(1) * (lambda(1) - lambdastart);
    end
    if lambdaend > lambda(n - 1)
        sum = sum + vals(n - 1) * (lambdaend - lambda(n - 1));
    end

    i = 1;
    while lambdastart > lambda(i+1) && i+1<n
        i = i+1;
    end
    while i + 1 < n && lambdaend >= lambda(i)
        seglambdastart = max(lambdastart, lambda(i));
        seglambdaend = min(lambdaend, lambda(i+1));
        sum = sum + 0.5 * (interp(lambda, seglambdastart, i, vals) + interp(lambda, seglambdaend, i, vals)) * (seglambdaend-seglambdastart);
        i = i + 1;
    end
    average = sum / (lambdaend - lambdastart);
end

function result =  interp(lambda, w, i, vals)
    result = lerp((w-lambda(i))/(lambda(i+1)-lambda(i)), vals(i), vals(i+1));
end

function val =  lerp(t, v0, v1)
    val = (1 - t) * v0 + t * v1;
end
