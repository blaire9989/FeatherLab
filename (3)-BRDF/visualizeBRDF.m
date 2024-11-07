function visualizeBRDF(type, name, dim, nExample, resolution, D65)

    % Load the noise functions for computing BRDF instances
    id = fopen("data" + num2str(type) + "/" + name + "/render/noise.binary", 'r');
    A = fread(id, [1 1], 'float');
    noise = fread(id, [1024, 10000], 'float');
    fclose(id);

    % Compute the smooth and 1-instance BRDF lobes for a few example incident directions
    brdf = cell(nExample + 1, 12);
    for d = 1 : 12
        x0 = 0;
        y0 = -0.88 + 0.08 * (d - 1);
        theta_i = asin(y0);
        phi_i = asin(x0 / cos(theta_i));
        [m1, p1, f1, m2, p2, f2, r, l] = evaluateParameters(type, name, dim, theta_i, phi_i);
        [smooth, signal, brdf{1, d}] = smoothBRDF(theta_i, phi_i, m1, p1, f1, m2, p2, f2, resolution, D65);
        for n = 1 : nExample
            brdf{n + 1, d} = individualBRDF(theta_i, phi_i, r, l, smooth, signal, A, noise(:, 2 * n - 1), noise(:, 2 * n), resolution, D65);
        end
    end
    
    % Output images
    mkdir("data" + num2str(type) + "/" + name + "/visual");
    for n = 0 : nExample
        collage = [brdf{n + 1, 1} brdf{n + 1, 5} brdf{n + 1, 9};
                   brdf{n + 1, 2} brdf{n + 1, 6} brdf{n + 1, 10};
                   brdf{n + 1, 3} brdf{n + 1, 7} brdf{n + 1, 11};
                   brdf{n + 1, 4} brdf{n + 1, 8} brdf{n + 1, 12}];
        imwrite(collage, "data" + num2str(type) + "/" + name + "/visual/brdf" + num2str(n) + ".jpg");
    end

end

function [m1, p1, f1, m2, p2, f2, r, l] = evaluateParameters(type, name, dim, theta_i, phi_i)

    % Load the parameter tables for bicubic interpolation
    id = fopen("data" + num2str(type) + "/" + name + "/render/parameters.binary", 'r');
    parameters = fread(id, [192 50 * dim * dim], 'float');
    fclose(id);
    m1 = zeros(50, 1);
    p1 = zeros(50, 1);
    f1 = zeros(50, 1);
    m2 = zeros(50, 1);
    p2 = zeros(50, 1);
    f2 = zeros(50, 1);
    r = zeros(50, 1);
    l = zeros(50, 1);

    % Evaluate relevant BRDF parameters for each wavelength
    for w = 1 : 50
        index1 = dim * dim * (w - 1) + 1;
        index2 = dim * dim * w;
        m1(w) = interpolate(parameters(1 : 16, index1 : index2), dim, cos(theta_i), cos(phi_i));
        m1(w) = max(m1(w), 0.001);
        m1(w) = min(m1(w), 1.000);
        p1(w) = interpolate(parameters(17 : 32, index1 : index2), dim, cos(theta_i), cos(phi_i));
        p1(w) = max(p1(w), -0.5 * pi);
        p1(w) = min(p1(w), 0.5 * pi);
        f1(w) = interpolate(parameters(33 : 48, index1 : index2), dim, cos(theta_i), cos(phi_i)) / (cos(theta_i) * cos(phi_i));
        m2(w) = interpolate(parameters(81 : 96, index1 : index2), dim, cos(theta_i), cos(phi_i));
        m2(w) = max(m2(w), 0.001);
        m2(w) = min(m2(w), 1.000);
        p2(w) = interpolate(parameters(97 : 112, index1 : index2), dim, cos(theta_i), cos(phi_i));
        p2(w) = max(p2(w), -0.5 * pi);
        p2(w) = min(p2(w), 0.5 * pi);
        f2(w) = interpolate(parameters(113 : 128, index1 : index2), dim, cos(theta_i), cos(phi_i)) / (cos(theta_i) * cos(phi_i));
        r(w) = interpolate(parameters(161 : 176, index1 : index2), dim, cos(theta_i), cos(phi_i));
        l(w) = interpolate(parameters(177 : 192, index1 : index2), dim, cos(theta_i), cos(phi_i));
    end

end

function q = interpolate(Q, dim, cos_theta, cos_phi)

    % Perform bicubic interpolation on prepared data
    x0 = cos_theta * dim;
    i = floor(x0) + 1;
    if i > dim
        i = dim;
    end
    s = x0 - (i - 1);
    y0 = cos_phi * dim;
    j = floor(y0) + 1;
    if j > dim
        j = dim;
    end
    t = y0 - (j - 1);
    index = dim * (i - 1) + j;
    interp = reshape(Q(:, index), [4, 4]);
    q = [1 s s^2 s^3] * interp * [1; t; t^2; t^3];

end

function [smooth, signal, RGBLobe] = smoothBRDF(theta_i, phi_i, m1, p1, f1, m2, p2, f2, resolution, D65)

    % Determine the incident direction
    flip1 = theta_i < 0;
    theta_i = abs(theta_i);
    flip2 = phi_i < 0;
    sigma = 0.04;
    smooth = zeros(resolution * resolution, 50);
    signal = zeros(resolution * resolution, 50);

    % Compute the smooth BRDF and statistics function (signal) for each wavelength
    for w = 1 : 50
        wi1 = [cos(theta_i) * sin(p1(w)); sin(theta_i); cos(theta_i) * cos(p1(w))];
        wi2 = [cos(theta_i) * sin(p2(w)); sin(theta_i); cos(theta_i) * cos(p2(w))];
        lobe1 = zeros(resolution, resolution);
        lobe2 = zeros(resolution, resolution);
        Gi1 = 2 / (1 + sqrt(1 + m1(w)^2 * (1 - wi1(3)^2) / wi1(3)^2));
        Gi2 = 2 / (1 + sqrt(1 + m2(w)^2 * (1 - wi2(3)^2) / wi2(3)^2));
        for i = 1 : resolution
            for j = 1 : resolution
                x0 = -1.0 + 1.0 / resolution + (i - 1) * 2.0 / resolution;
                y0 = -1.0 + 1.0 / resolution + (j - 1) * 2.0 / resolution;
                if x0 * x0 + y0 * y0 > 1
                    continue;
                end
                wo = [x0; y0; sqrt(1 - x0 * x0 - y0 * y0)];
                Go1 = 2 / (1 + sqrt(1 + m1(w)^2 * (1 - wo(3)^2) / wo(3)^2));
                Go2 = 2 / (1 + sqrt(1 + m2(w)^2 * (1 - wo(3)^2) / wo(3)^2));
                h1 = (wi1 + wo) / norm(wi1 + wo);
                h2 = (wi2 + wo) / norm(wi2 + wo);
                D11 = 1 / pi * m1(w)^2 / (m1(w)^2 + (1 - m1(w)^2) * h1(1)^2)^2;
                D21 = exp(-h1(2)^2 / sigma^2);
                lobe1(resolution + 1 - j, i) = f1(w) * D11 * D21 * Gi1 * Go1 / (4 * wi1(3));
                D12 = 1 / pi * m2(w)^2 / (m2(w)^2 + (1 - m2(w)^2) * h2(1)^2)^2;
                D22 = exp(-h2(2)^2 / sigma^2);
                lobe2(resolution + 1 - j, i) = f2(w) * D12 * D22 * Gi2 * Go2 / (4 * wi2(3));
                if abs(lobe1(resolution + 1 - j, i)) < 1e-20
                    lobe1(resolution + 1 - j, i) = 0;
                end
                if abs(lobe2(resolution + 1 - j, i)) < 1e-20
                    lobe2(resolution + 1 - j, i) = 0;
                end
            end
        end
        if flip1
            lobe1 = flip(lobe1, 1);
            lobe2 = flip(lobe2, 1);
        end
        if flip2
            lobe1 = flip(lobe1, 2);
            lobe2 = flip(lobe2, 2);
        end
        smooth(:, w) = reshape(lobe1, [resolution * resolution, 1]);
        signal(:, w) = reshape(lobe2, [resolution * resolution, 1]);
    end

    % Convert the smooth BRDF in the RGB (gamma-corrected)
    spec = ToXYZ(smooth, 50, 400, 700, D65);
    spec = XYZToRGB(spec);
    RGBLobe = zeros(resolution, resolution, 3);
    RGBLobe(:, :, 1) = reshape(spec(:, 1), [resolution resolution]);
    RGBLobe(:, :, 2) = reshape(spec(:, 2), [resolution resolution]);
    RGBLobe(:, :, 3) = reshape(spec(:, 3), [resolution resolution]);
    RGBLobe(RGBLobe < 0) = 0;
    RGBLobe = RGBLobe .^ (1 / 2.2);
    for i = 1 : resolution
        for j = 1 : resolution
            x0 = -1.0 + 1.0 / resolution + (i - 1) * 2.0 / resolution;
            y0 = -1.0 + 1.0 / resolution + (j - 1) * 2.0 / resolution;
            if x0 * x0 + y0 * y0 > 1
                RGBLobe(resolution + 1 - j, i, 1) = 0.8;
                RGBLobe(resolution + 1 - j, i, 2) = 0.8;
                RGBLobe(resolution + 1 - j, i, 3) = 0.8;
            end
        end
    end

end

function RGBLobe = individualBRDF(theta_i, phi_i, r, l, smooth, signal, A, noiseR, noiseI, resolution, D65)

    % Some constants
    c0 = 299792458.0;
    mu = 4e-7 * pi;
    eps = 1 / (mu * c0 * c0);
    scale = 0.5 * sqrt(eps / mu);
    instance = zeros(resolution * resolution, 50);

    % Compute the BRDF instance at each wavelength
    for w = 1 : 50
        smoothW = reshape(smooth(:, w), [resolution, resolution]);
        signalW = reshape(signal(:, w), [resolution, resolution]);
        lobe = zeros(resolution, resolution);
        valueR = zeros(resolution, 1);
        valueI = zeros(resolution, 1);
        for i = 1 : resolution
            x0 = -1.0 + 1.0 / resolution + (i - 1) * 2.0 / resolution;
            valueR(i) = evaluateNoise(noiseR, l(w), min(l), cos(theta_i) * sin(phi_i), x0, A);
            valueI(i) = evaluateNoise(noiseI, l(w), min(l), cos(theta_i) * sin(phi_i), x0, A);
        end
        boundL1 = 1;
        while boundL1 <= resolution && valueR(boundL1) == -Inf
            boundL1 = boundL1 + 1;
        end
        boundR1 = resolution;
        while boundR1 >= 1 && valueR(boundR1) == Inf
            boundR1 = boundR1 - 1;
        end
        for j = 1 : resolution
            boundL2 = resolution + 1;
            boundR2 = 0;
            for i = 1 : resolution
                x0 = -1.0 + 1.0 / resolution + (i - 1) * 2.0 / resolution;
                y0 = -1.0 + 1.0 / resolution + (j - 1) * 2.0 / resolution;
                if x0 * x0 + y0 * y0 > 1
                    continue;
                end
                boundL2 = min(i, boundL2);
                boundR2 = max(i, boundR2);
                if i < boundL1 || i > boundR1
                    lobe(resolution + 1 - j, i) = smoothW(resolution + 1 - j, i);
                    continue;
                end
                [muR1, muI1, sdR1, sdI1, muR2, muI2, sdR2, sdI2] = singlePointStatistics(i, j, r(w), scale, smoothW, signalW, resolution);
                fieldR1 = sdR1 * valueR(i) + muR1;
                fieldI1 = sdI1 * valueI(i) + muI1;
                fieldR2 = -sdR2 * valueR(i) + muR2;
                fieldI2 = -sdI2 * valueI(i) + muI2;
                lobe(resolution + 1 - j, i) = 0.5 * scale * (fieldR1^2 + fieldI1^2 + fieldR2^2 + fieldI2^2);
            end
            lobe(resolution + 1 - j, :) = linearBlend(lobe(resolution + 1 - j, :), boundL1, boundR1, boundL2, boundR2, 20);
        end
        instance(:, w) = reshape(lobe, [resolution * resolution, 1]);
    end

    % Convert the individual BRDF in the RGB (gamma-corrected)
    spec = ToXYZ(instance, 50, 400, 700, D65);
    spec = XYZToRGB(spec);
    RGBLobe = zeros(resolution, resolution, 3);
    RGBLobe(:, :, 1) = reshape(spec(:, 1), [resolution resolution]);
    RGBLobe(:, :, 2) = reshape(spec(:, 2), [resolution resolution]);
    RGBLobe(:, :, 3) = reshape(spec(:, 3), [resolution resolution]);
    RGBLobe(RGBLobe < 0) = 0;
    RGBLobe = RGBLobe .^ (1 / 2.2);
    for i = 1 : resolution
        for j = 1 : resolution
            x0 = -1.0 + 1.0 / resolution + (i - 1) * 2.0 / resolution;
            y0 = -1.0 + 1.0 / resolution + (j - 1) * 2.0 / resolution;
            if x0 * x0 + y0 * y0 > 1
                RGBLobe(resolution + 1 - j, i, 1) = 0.8;
                RGBLobe(resolution + 1 - j, i, 2) = 0.8;
                RGBLobe(resolution + 1 - j, i, 3) = 0.8;
            end
        end
    end

end

function value = evaluateNoise(noise, l, lmin, xi, xo, A)

    % Scale and translate the domain of the given noise function
    N = length(noise);
    xvec = l * linspace(-1.0 - 1.0 / N, 1.0 + 1.0 / N, N + 2);
    noise = [0; noise; 0];
    x0 = xi + xo;

    % Interpolate the noise function at the queried point, preserving the standard deviation
    if x0 < -lmin
        value = -Inf;
    elseif x0 > lmin
        value = Inf;
    else
        index = 1;
        while index <= length(xvec) && xvec(index) <= x0
            index = index + 1;
        end
        s = (x0 - xvec(index - 1)) / (xvec(index) - xvec(index - 1));
        r1 = (1 - s) / sqrt((1 - s)^2 + s^2 + 2 * (1 - s) * s * A);
        r2 = s / sqrt((1 - s)^2 + s^2 + 2 * (1 - s) * s * A);
        value = r1 * noise(index - 1) + r2 * noise(index);
    end

end

function [muR1, muI1, sdR1, sdI1, muR2, muI2, sdR2, sdI2] = singlePointStatistics(i, j, r, scale, smooth, signal, resolution)

    % Compute the mean and standard deviation of the scattered field values
    intensity1 = 2.0 * r * smooth(resolution + 1 - j, i);
    signal1 = 2.0 * r * signal(resolution + 1 - j, i);
    signal1 = min(intensity1, signal1);
    var1 = intensity1 - signal1;
    intensity2 = 2.0 * (1 - r) * smooth(resolution + 1 - j, i);
    signal2 = 2.0 * (1 - r) * signal(resolution + 1 - j, i);
    signal2 = min(intensity2, signal2);
    var2 = intensity2 - signal2;
    muR1 = sqrt(signal1 / scale);
    muI1 = 0;
    sdR1 = sqrt(0.5 * var1 / scale);
    sdI1 = sqrt(0.5 * var1 / scale);
    muR2 = sqrt(signal2 / scale);
    muI2 = 0;
    sdR2 = sqrt(0.5 * var2 / scale);
    sdI2 = sqrt(0.5 * var2 / scale);

end

function slice = linearBlend(slice, boundL1, boundR1, boundL2, boundR2, grace)

    % Blend between the BRDF slices with and without adding in noise (for visualization only)
    if boundL2 < boundL1
        index = max(boundL1 - grace, boundL2);
        for i = index : boundL1 - 1
            s = abs(i - boundL1) / grace;
            slice(i) = (1 - s) * slice(boundL1) + s * slice(i);
        end
    end
    if boundR2 > boundR1
        index = min(boundR1 + grace, boundR2);
        for i = boundR1 + 1 : index
            s = abs(i - boundR1) / grace;
            slice(i) = (1 - s) * slice(boundR1) + s * slice(i);
        end
    end

end