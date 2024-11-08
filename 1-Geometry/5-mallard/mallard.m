function [xyvals, info, regularity] = mallard(width, thickness, d)
addpath('../utils/');

    % Generate the outer contour and only model the barbule top surface
    N = round(width / d);
    if rem(N, 2) == 1
        N = N + 1;
    end
    xtop = transpose(linspace(-0.5 * width, 0.5 * width, N + 1));
    ytop = zeros(N + 1, 1);
    xbot = transpose(linspace(-0.5 * width, 0.5 * width, N + 1));
    ybot = -thickness * ones(N + 1, 1);
    contour = fullContour(xtop, ytop, xbot, ybot, d, 2.0 + 1.0 * rand);
    perturb = roughness(1.2, -width, width, -0.3, 0.3, size(contour, 1), false);
    contour(:, 2) = contour(:, 2) + perturb(1 : size(contour, 1));
    [xyvals, ~] = splitGeometry(contour, 0.2);
    xyvals = xyvals * width / (max(xyvals(:, 1)) - min(xyvals(:, 1)));
    xshift = 0.5 * (min(xyvals(:, 1)) + max(xyvals(:, 1)));
    yshift = mean(xyvals(:, 2));
    xyvals(:, 1) = xyvals(:, 1) - xshift;
    xyvals(:, 2) = xyvals(:, 2) - yshift;
    info = [0 size(xyvals, 1)];
    regularity = ones(size(xyvals, 1) - 2, 1);

end