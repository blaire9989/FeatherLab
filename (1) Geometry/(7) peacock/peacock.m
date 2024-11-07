function [xyvals, info, regularity] = peacock(width, d)
addpath('../utils/');

    % Generate the slightly curved barbule top surface
    slope = tan(pi / 12 + pi / 15 * rand);
    phi = pi / 6 + pi / 4 * rand;
    a0 = 0.5 * width / cos(phi);
    b0 = a0 * tan(phi) * slope;
    L = ellipseArc(a0, b0, phi, pi - phi);
    N = ceil(L / d);
    if rem(N, 2) == 1
        N = N + 1;
    end
    d0 = L / N;
    xtop = zeros(N + 1, 1);
    ytop = zeros(N + 1, 1);
    count = 1;
    phi0 = pi - phi;
    while count <= N + 1
        xtop(count) = a0 * cos(phi0);
        ytop(count) = b0 * sin(phi0);
        dPhi = d0 / sqrt(a0^2 * sin(phi0)^2 + b0^2 * cos(phi0)^2);
        phi0 = phi0 - dPhi;
        count = count + 1;
    end
    perturb = roughness(1.2, -0.5 * width, 0.5 * width, -0.3, 0.3, N + 1, true);
    ytop = ytop + perturb(1 : N + 1);
    ytop = ytop - mean(ytop);
    xyvals = [xtop ytop];
    info = [0 size(xyvals, 1)];
    regularity = ones(size(xyvals, 1) - 2, 1);

end