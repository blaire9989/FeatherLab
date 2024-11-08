function constructMaterial(type, name, dim)

    % Obtain the reference autocorrelation length
    id = fopen("data" + num2str(type) + "/" + name + "/render/param_400_400.binary", 'r');
    param = fread(id, [12 1], 'double');
    fclose(id);
    l0 = param(12);

    % Extract the parameters for each wavelength and compute bicubic interpolation tables
    M1 = processM(type, name, dim, 1);
    P1 = processP(type, name, dim, 2);
    F1_s = processF(type, name, dim, 3);
    F1_m = processF(type, name, dim, 4);
    F1_l = processF(type, name, dim, 5);
    M2 = processM(type, name, dim, 6);
    P2 = processP(type, name, dim, 7);
    F2_s = processF(type, name, dim, 8);
    F2_m = processF(type, name, dim, 9);
    F2_l = processF(type, name, dim, 10);
    R = processR(type, name, dim, 11);
    L = processL(type, name, dim, 12, l0);
    parameters = [M1; P1; F1_s; F1_m; F1_l; M2; P2; F2_s; F2_m; F2_l; R; L];
    for wvl = 400 : 6 : 694
        for index = 1 : 400
            delete("data" + num2str(type) + "/" + name + "/render/param_" + num2str(wvl) + "_" + num2str(index) + ".binary");
        end
    end
    id = fopen("data" + num2str(type) + "/" + name + "/render/parameters.binary", 'w');
    fwrite(id, single(parameters), 'float');
    fclose(id);

    % Generate noise tiles for individual variations
    resolution = 1024;
    nInstance = 50;
    id = fopen("data" + num2str(type) + "/" + name + "/render/noise.binary", 'r');
    noise = fread(id, [resolution 4 * nInstance], 'float');
    fclose(id);
    acv = computeAverageACV(noise, resolution, nInstance);
    A = acv(2);
    psd = abs(fft(acv) / resolution);
    noise = generateNoise(psd, resolution, 10000);
    delete("data" + num2str(type) + "/" + name + "/render/noise.binary");
    id = fopen("data" + num2str(type) + "/" + name + "/render/noise.binary", 'w');
    fwrite(id, single(A), 'float');
    fwrite(id, single(noise), 'float');
    fclose('all');

end

function M = processM(type, name, dim, pos)

    % Extract the roughness m parameter for each simulated wavelength and build bicubic interpolation tables
    M = zeros(16, 50 * dim * dim);
    for i = 1 : 50
        wvl = 400 + 6 * (i - 1);
        m_raw = zeros(dim + 1, dim + 1);
        for j = 1 : dim
            for k = 1 : dim
                index = dim * (j - 1) + k;
                id = fopen("data" + num2str(type) + "/" + name + "/render/param_" + num2str(wvl) + "_" + num2str(index) + ".binary", 'r');
                param = fread(id, [12 1], 'double');
                m_raw(j + 1, k + 1) = param(pos);
                fclose(id);
            end
        end
        for j = 1 : dim
            value = 2 * m_raw(j + 1, 2) - m_raw(j + 1, 3);
            value = max(0.0, value);
            value = min(1.0, value);
            m_raw(j + 1, 1) = value;
        end
        m_raw(1, :) = 2 * mean(m_raw(2, :)) - mean(m_raw(3, :));
        m = zeros(dim + 1, dim + 1);
        for mx = 1 : dim + 1
            for my = 1 : dim + 1
                value = 0;
                weight = 0;
                for nx = 1 : dim + 1
                    for ny = 1 : dim + 1
                        if abs(nx - mx) >= 3 && abs(ny - my) >= 3
                            continue;
                        end
                        g = exp(-(nx - mx)^2) * exp(-(ny - my)^2);
                        value = value + g * m_raw(nx, ny);
                        weight = weight + g;
                    end
                end
                m(mx, my) = value / weight;
            end
        end
        M(:, dim * dim * (i - 1) + 1 : dim * dim * i) = buildBicubicTable(m, dim);
    end

end

function P = processP(type, name, dim, pos)

    % Extract the effective phi_i parameter for each simulated wavelength and build bicubic interpolation tables
    P = zeros(16, 50 * dim * dim);
    for i = 1 : 50
        wvl = 400 + 6 * (i - 1);
        p = zeros(dim + 1, dim + 1);
        for j = 1 : dim
            for k = 1 : dim
                index = dim * (j - 1) + k;
                id = fopen("data" + num2str(type) + "/" + name + "/render/param_" + num2str(wvl) + "_" + num2str(index) + ".binary", 'r');
                param = fread(id, [12 1], 'double');
                p(j + 1, k + 1) = param(pos);
                fclose(id);
            end
        end
        p(:, 1) = pi / 2;
        for j = 1 : dim
            p(1, j + 1) = acos(j / dim);
        end
        P(:, dim * dim * (i - 1) + 1 : dim * dim * i) = buildBicubicTable(p, dim);
    end

end

function F = processF(type, name, dim, pos)

    % Extract the effective Fresnel parameter for each simulated wavelength and build bicubic interpolation tables
    F = zeros(16, 50 * dim * dim);
    for i = 1 : 50
        wvl = 400 + 6 * (i - 1);
        f = zeros(dim + 1, dim + 1);
        for j = 1 : dim
            for k = 1 : dim
                index = dim * (j - 1) + k;
                id = fopen("data" + num2str(type) + "/" + name + "/render/param_" + num2str(wvl) + "_" + num2str(index) + ".binary", 'r');
                param = fread(id, [12 1], 'double');
                f(j + 1, k + 1) = param(pos);
                fclose(id);
            end
        end
        F(:, dim * dim * (i - 1) + 1 : dim * dim * i) = buildBicubicTable(f, dim);
    end

end

function R = processR(type, name, dim, pos)

    % Extract the polarization ratio parameter for each simulated wavelength and build bicubic interpolation tables
    R = zeros(16, 50 * dim * dim);
    for i = 1 : 50
        wvl = 400 + 6 * (i - 1);
        r = zeros(dim + 1, dim + 1);
        for j = 1 : dim
            for k = 1 : dim
                index = dim * (j - 1) + k;
                id = fopen("data" + num2str(type) + "/" + name + "/render/param_" + num2str(wvl) + "_" + num2str(index) + ".binary", 'r');
                param = fread(id, [12 1], 'double');
                r(j + 1, k + 1) = param(pos);
                fclose(id);
            end
        end
        for j = 1 : dim
            r(j + 1, 1) = 2 * r(j + 1, 2) - r(j + 1, 3);
        end
        for j = 1 : dim + 1
            r(1, j) = 2 * r(2, j) - r(3, j);
        end
        R(:, dim * dim * (i - 1) + 1 : dim * dim * i) = buildBicubicTable(r, dim);
    end

end

function L = processL(type, name, dim, pos, l0)

    % Extract the roughness m parameter for each simulated wavelength and build bicubic interpolation tables
    L = zeros(16, 50 * dim * dim);
    for i = 1 : 50
        wvl = 400 + 6 * (i - 1);
        l_raw = zeros(dim + 1, dim + 1);
        for j = 1 : dim
            for k = 1 : dim
                index = dim * (j - 1) + k;
                id = fopen("data" + num2str(type) + "/" + name + "/render/param_" + num2str(wvl) + "_" + num2str(index) + ".binary", 'r');
                param = fread(id, [12 1], 'double');
                l_raw(j + 1, k + 1) = param(pos) / l0;
                fclose(id);
            end
        end
        for j = 1 : dim
            l_raw(j + 1, 1) = 2 * l_raw(j + 1, 2) - l_raw(j + 1, 3);
        end
        l_raw(1, :) = 2 * mean(l_raw(2, :)) - mean(l_raw(3, :));
        l = zeros(dim + 1, dim + 1);
        for mx = 1 : dim + 1
            for my = 1 : dim + 1
                value = 0;
                weight = 0;
                for nx = 1 : dim + 1
                    for ny = 1 : dim + 1
                        if abs(nx - mx) >= 3 && abs(ny - my) >= 3
                            continue;
                        end
                        g = exp(-(nx - mx)^2) * exp(-(ny - my)^2);
                        value = value + g * l_raw(nx, ny);
                        weight = weight + g;
                    end
                end
                l(mx, my) = value / weight;
            end
        end
        L(:, dim * dim * (i - 1) + 1 : dim * dim * i) = buildBicubicTable(l, dim);
    end

end

function Q = buildBicubicTable(q, dim)

    % Compute function values and derivatives at each corner point
    qx = zeros(dim + 1, dim + 1);
    qy = zeros(dim + 1, dim + 1);
    qxy = zeros(dim + 1, dim + 1);
    for i = 1 : dim + 1
        for j = 1 : dim + 1
            i1 = max(1, i - 1);
            i2 = min(dim + 1, i + 1);
            j1 = max(1, j - 1);
            j2 = min(dim + 1, j + 1);
            qx(i, j) = (q(i2, j) - q(i1, j)) / (i2 - i1);
            qy(i, j) = (q(i, j2) - q(i, j1)) / (j2 - j1);
            qxy(i, j) = (q(i1, j1) + q(i2, j2) - q(i1, j2) - q(i2, j1)) / ((i2 - i1) * (j2 - j1));
        end
    end

    % Compute the bicubic coefficients a00-a33 for the dim * dim blocks
    Q = zeros(16, dim * dim);
    C = [1 0 0 0; 0 0 1 0; -3 3 -2 -1; 2 -2 1 1];
    for i = 1 : dim
        for j = 1 : dim
            A = [q(i, j)      q(i, j + 1)      qy(i, j)      qy(i, j + 1);
                 q(i + 1, j)  q(i + 1, j + 1)  qy(i + 1, j)  qy(i + 1, j + 1);
                 qx(i, j)     qx(i, j + 1)     qxy(i, j)     qxy(i, j + 1);
                 qx(i + 1, j) qx(i + 1, j + 1) qxy(i + 1, j) qxy(i + 1, j + 1)];
            B = C * A * transpose(C);
            index = dim * (i - 1) + j;
            Q(:, index) = reshape(B, [16, 1]);
        end
    end

end

function acv = computeAverageACV(noise, resolution, nInstance)

    % The noise function is periodic since it evaluates to 0 at both ends
    acv = zeros(resolution / 2, 4 * nInstance);
    for i = 1 : 4 * nInstance
        average = mean(noise(:, i));
        for j = 0 : resolution / 2 - 1
            total = 0;
            for k = 1 : resolution
                index = k + j;
                if index > resolution
                    index = index - resolution;
                end
                total = total + (noise(k, i) - average) * (noise(index, i) - average);
            end
            acv(j + 1, i) = total / (resolution - 1);
        end
    end

    % Compute the average autocorrelation function
    acv = mean(acv, 2);
    acv = [acv; 0; flip(acv(2 : end))];
    acv(resolution / 2 + 1) = 0.5 * acv(resolution / 2) + 0.5 * acv(resolution / 2 + 2);
    acv = acv / acv(1);

end

function noise = generateNoise(psd, resolution, N)

    % Generate a tile of many noise functions
    phase = normrnd(0, sqrt(0.5), resolution / 2 - 1, 2 * N);
    noise = zeros(resolution, N);
    for i = 1 : N
        phaseR = phase(:, 2 * i - 1);
        phaseI = phase(:, 2 * i);
        spec = zeros(resolution, 1);
        spec(1) = sqrt(psd(1));
        spec(2 : resolution / 2) = sqrt(psd(2 : resolution / 2)) .* (phaseR + 1j * phaseI);
        spec(resolution / 2 + 2 : end) = conj(flip(spec(2 : resolution / 2)));
        noise(:, i) = real(resolution * ifft(spec));
    end

end