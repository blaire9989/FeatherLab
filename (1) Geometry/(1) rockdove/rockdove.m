function [xyvals, info, mRatio] = rockdove(width, thickness, cortex, d)
addpath('../utils/');

    % Generate model layer 1: top outer cortex
    slope = tan(pi / 10 + pi / 15 * rand);
    phi = pi / 6 + pi / 4 * rand;
    a0 = 0.5 * width / cos(phi);
    b0 = a0 * tan(phi) * slope;
    L = ellipseArc(a0, b0, phi, pi - phi);
    N = ceil(L / d);
    if rem(N, 2) == 1
        N = N + 1;
    end
    d0 = L / N;
    xlayer1 = zeros(N + 1, 1);
    ylayer1 = zeros(N + 1, 1);
    count = 1;
    phi0 = pi - phi;
    while count <= N + 1
        xlayer1(count) = a0 * cos(phi0);
        ylayer1(count) = b0 * sin(phi0);
        dPhi = d0 / sqrt(a0^2 * sin(phi0)^2 + b0^2 * cos(phi0)^2);
        phi0 = phi0 - dPhi;
        count = count + 1;
    end
    perturb = roughness(1.0, -0.5 * width, 0.5 * width, -0.08, 0.08, length(ylayer1), false);
    ylayer1 = ylayer1 + perturb(1 : length(ylayer1));
    ylayer1 = ylayer1 - mean(ylayer1);

    % Generate model layer 2: top inner cortex
    xlayer2 = zeros(length(xlayer1), 1);
    ylayer2 = zeros(length(ylayer1), 1);
    for i = 1 : length(xlayer2)
        if i == 1
            tvec = [xlayer1(2) - xlayer1(1) ylayer1(2) - ylayer1(1)];
        elseif i == length(xlayer2)
            tvec = [xlayer1(end) - xlayer1(end - 1) ylayer1(end) - ylayer1(end - 1)];
        else
            tvec = [xlayer1(i + 1) - xlayer1(i - 1) ylayer1(i + 1) - ylayer1(i - 1)];
        end
        tvec = tvec / norm(tvec);
        xlayer2(i) = xlayer1(i) + cortex * tvec(2);
        ylayer2(i) = ylayer1(i) - cortex * tvec(1);
    end

    % Generate model layer 3: fictitious boundary between air and melanin in a rock dove barbule
    airlayer = roughness(0.2, -0.5 * width, 0.5 * width, 0.25, 0.75, length(ylayer2), false);
    airlayer = airlayer(1 : length(ylayer2));
    xlayer3 = zeros(length(xlayer2), 1);
    ylayer3 = zeros(length(ylayer2), 1);
    for i = 1 : length(xlayer3)
        if i == 1
            tvec = [xlayer2(2) - xlayer2(1) ylayer2(2) - ylayer2(1)];
        elseif i == length(xlayer3)
            tvec = [xlayer2(end) - xlayer2(end - 1) ylayer2(end) - ylayer2(end - 1)];
        else
            tvec = [xlayer2(i + 1) - xlayer2(i - 1) ylayer2(i + 1) - ylayer2(i - 1)];
        end
        tvec = tvec / norm(tvec);
        xlayer3(i) = xlayer2(i) + airlayer(i) * tvec(2);
        ylayer3(i) = ylayer2(i) - airlayer(i) * tvec(1);
    end

    % Generate model layer 4: lower inner cortex
    slope = 0.9 * slope;
    phi = pi / 6 + pi / 4 * rand;
    a0 = 0.5 * width / cos(phi);
    b0 = a0 * tan(phi) * slope;
    L = ellipseArc(a0, b0, phi, pi - phi);
    N = ceil(L / d);
    if rem(N, 2) == 1
        N = N + 1;
    end
    d0 = L / N;
    xlayer4 = zeros(N + 1, 1);
    ylayer4 = zeros(N + 1, 1);
    count = 1;
    phi0 = pi - phi;
    while count <= N + 1
        xlayer4(count) = a0 * cos(phi0);
        ylayer4(count) = b0 * sin(phi0);
        dPhi = d0 / sqrt(a0^2 * sin(phi0)^2 + b0^2 * cos(phi0)^2);
        phi0 = phi0 - dPhi;
        count = count + 1;
    end
    perturb = roughness(0.7, -0.45 * width, 0.45 * width, -0.07, 0.07, length(ylayer4), false);
    ylayer4 = ylayer4 + perturb(1 : length(ylayer4));
    ylayer4 = ylayer4 - mean(ylayer4) - thickness;

    % Generate model layer 5: lower outer cortex
    xlayer5 = zeros(length(xlayer4), 1);
    ylayer5 = zeros(length(ylayer4), 1);
    for i = 1 : length(xlayer5)
        if i == 1
            tvec = [xlayer4(2) - xlayer4(1) ylayer4(2) - ylayer4(1)];
        elseif i == length(xlayer5)
            tvec = [xlayer4(end) - xlayer4(end - 1) ylayer4(end) - ylayer4(end - 1)];
        else
            tvec = [xlayer4(i + 1) - xlayer4(i - 1) ylayer4(i + 1) - ylayer4(i - 1)];
        end
        tvec = tvec / norm(tvec);
        xlayer5(i) = xlayer4(i) + cortex * tvec(2);
        ylayer5(i) = ylayer4(i) - cortex * tvec(1);
    end

    % Put together the geometric model
    info = [0 size(xlayer1, 1)];
    xyvals = [xlayer1 ylayer1];
    info = [info; size(xyvals, 1) size(xlayer2, 1)];
    xyvals = [xyvals; xlayer2 ylayer2];
    info = [info; size(xyvals, 1) size(xlayer3, 1)];
    xyvals = [xyvals; xlayer3 ylayer3];
    info = [info; size(xyvals, 1) size(xlayer4, 1)];
    xyvals = [xyvals; xlayer4 ylayer4];
    info = [info; size(xyvals, 1) size(xlayer5, 1)];
    xyvals = [xyvals; xlayer5 ylayer5];

    % Model the percentage of melanin in the interior
    total = size(xlayer1, 1) - 2;
    if rem(total, 2) == 1
        total = total + 1;
    end
    mRatio = roughness(1.0, -0.5 * width, 0.5 * width, 0.8, 1.0, total, false);
    mRatio = mRatio(1 : size(xlayer1, 1) - 2);

end