function [xyvals, info, mRatio] = starling(width, thickness, cortex, Dm, d)
addpath('../utils/');

    % Generate the top and bottom surfaces
    N = round(width / d);
    if rem(N, 2) == 1
        N = N + 1;
    end
    xtop = transpose(linspace(-0.5 * width, 0.5 * width, N + 1));
    ytop = zeros(N + 1, 1);
    xbot = transpose(linspace(-0.5 * width, 0.5 * width, N + 1));
    ybot = -thickness * ones(N + 1, 1);
    contour = fullContour(xtop, ytop, xbot, ybot, d, 2.0 + 1.0 * rand);
    perturb = roughness(1.5, -width, width, -0.3, 0.3, size(contour, 1), false);
    contour(:, 2) = contour(:, 2) + perturb(1 : size(contour, 1));
    [top, bot] = splitGeometry(contour, 0.2);
    bot = bot * width / (max(top(:, 1)) - min(top(:, 1)));
    top = top * width / (max(top(:, 1)) - min(top(:, 1)));
    xshift = 0.5 * (min(top(:, 1)) + max(top(:, 1)));
    yshift = mean(top(:, 2));
    top(:, 1) = top(:, 1) - xshift;
    top(:, 2) = top(:, 2) - yshift;
    bot(:, 1) = bot(:, 1) - xshift;
    bot(:, 2) = bot(:, 2) - mean(bot(:, 2)) - thickness;
   
    % Locate melanosomes near the barbule top surface
    guide1 = zeros(size(top, 1), 2);
    for i = 1 : size(top, 1)
        if i == 1
            tvec = [top(2, 1) - top(1, 1) top(2, 2) - top(1, 2)];
        elseif i == size(top, 1)
            tvec = [top(end, 1) - top(end - 1, 1) top(end, 2) - top(end - 1, 2)];
        else
            tvec = [top(i + 1, 1) - top(i - 1, 1) top(i + 1, 2) - top(i - 1, 2)];
        end
        tvec = tvec / norm(tvec);
        guide1(i, :) = top(i, :) - (cortex + 0.5 * Dm) * [-tvec(2) tvec(1)];
    end
    slopeL = (guide1(2, 2) - guide1(1, 2)) / (guide1(2, 1) - guide1(1, 1));
    x0 = guide1(1, 1) - 2.0;
    y0 = guide1(1, 2) - 2.0 * slopeL;
    n = ceil(2.0 * sqrt(1 + slopeL * slopeL) / d);
    xExtra = linspace(x0, guide1(1, 1), n + 1);
    yExtra = linspace(y0, guide1(1, 2), n + 1);
    guide1 = [transpose(xExtra(1 : n)) transpose(yExtra(1 : n)); guide1];
    slopeR = (guide1(end, 2) - guide1(end - 1, 2)) / (guide1(end, 1) - guide1(end - 1, 1));
    x0 = guide1(end, 1) + 2.0;
    y0 = guide1(end, 2) + 2.0 * slopeR;
    n = ceil(2.0 * sqrt(1 + slopeR * slopeR) / d);
    xExtra = linspace(guide1(end, 1), x0, n + 1);
    yExtra = linspace(guide1(end, 2), y0, n + 1);
    guide1 = [guide1; transpose(xExtra(2 : n + 1)) transpose(yExtra(2 : n + 1))];

    % Locate melanosomes near the barbule bottom surface
    guide2 = zeros(size(bot, 1), 2);
    for i = 1 : size(bot, 1)
        if i == 1
            tvec = [bot(2, 1) - bot(1, 1) bot(2, 2) - bot(1, 2)];
        elseif i == size(bot, 1)
            tvec = [bot(end, 1) - bot(end - 1, 1) bot(end, 2) - bot(end - 1, 2)];
        else
            tvec = [bot(i + 1, 1) - bot(i - 1, 1) bot(i + 1, 2) - bot(i - 1, 2)];
        end
        tvec = tvec / norm(tvec);
        guide2(i, :) = bot(i, :) + (cortex + 0.5 * Dm) * [-tvec(2) tvec(1)];
    end
    slopeL = (guide2(2, 2) - guide2(1, 2)) / (guide2(2, 1) - guide2(1, 1));
    x0 = guide2(1, 1) - 2.0;
    y0 = guide2(1, 2) - 2.0 * slopeL;
    n = ceil(2.0 * sqrt(1 + slopeL * slopeL) / d);
    xExtra = linspace(x0, guide2(1, 1), n + 1);
    yExtra = linspace(y0, guide2(1, 2), n + 1);
    guide2 = [transpose(xExtra(1 : n)) transpose(yExtra(1 : n)); guide2];
    slopeR = (guide2(end, 2) - guide2(end - 1, 2)) / (guide2(end, 1) - guide2(end - 1, 1));
    x0 = guide2(end, 1) + 2.0;
    y0 = guide2(end, 2) + 2.0 * slopeR;
    n = ceil(2.0 * sqrt(1 + slopeR * slopeR) / d);
    xExtra = linspace(guide2(end, 1), x0, n + 1);
    yExtra = linspace(guide2(end, 2), y0, n + 1);
    guide2 = [guide2; transpose(xExtra(2 : n + 1)) transpose(yExtra(2 : n + 1))];

    % Generate layers of melanosomes
    [xc1, yc1, vD1] = growInterior(guide1, Dm, 0.2);
    [xc2, yc2, vD2] = growInterior(guide2, Dm, 0.2);
    melanin1 = generateLayer(xc1, yc1, vD1, d);
    melanin2 = generateLayer(xc2, yc2, vD2, d);
    [layer1, layer2] = splitGeometry(melanin1, 0.1);
    index1 = 1;
    while layer1(index1, 1) < -0.5 * width
        index1 = index1 + 1;
    end
    index2 = size(layer1, 1);
    while layer1(index2, 1) > 0.5 * width
        index2 = index2 - 1;
    end
    layer1 = layer1(index1 : index2, :);
    index1 = 1;
    while layer2(index1, 1) < -0.5 * width
        index1 = index1 + 1;
    end
    index2 = size(layer2, 1);
    while layer2(index2, 1) > 0.5 * width
        index2 = index2 - 1;
    end
    layer2 = layer2(index1 : index2, :);
    [layer3, layer4] = splitGeometry(melanin2, 0.1);
    index1 = 1;
    while layer3(index1, 1) < -0.5 * width
        index1 = index1 + 1;
    end
    index2 = size(layer3, 1);
    while layer3(index2, 1) > 0.5 * width
        index2 = index2 - 1;
    end
    layer3 = layer3(index1 : index2, :);
    index1 = 1;
    while layer4(index1, 1) < -0.5 * width
        index1 = index1 + 1;
    end
    index2 = size(layer4, 1);
    while layer4(index2, 1) > 0.5 * width
        index2 = index2 - 1;
    end
    layer4 = layer4(index1 : index2, :);

    % Put together the geometric model
    info = [0 size(top, 1)];
    xyvals = top;
    info = [info; size(xyvals, 1) size(layer1, 1)];
    xyvals = [xyvals; layer1];
    info = [info; size(xyvals, 1) size(layer2, 1)];
    xyvals = [xyvals; layer2];
    info = [info; size(xyvals, 1) size(layer3, 1)];
    xyvals = [xyvals; layer3];
    info = [info; size(xyvals, 1) size(layer4, 1)];
    xyvals = [xyvals; layer4];
    info = [info; size(xyvals, 1) size(bot, 1)];
    xyvals = [xyvals; bot];

    % Model the percentage of melanin in the interior
    total = size(top, 1) - 2;
    if rem(total, 2) == 1
        total = total + 1;
    end
    mRatio = roughness(1.0, -0.5 * width, 0.5 * width, 0.3, 0.7, total, false);
    mRatio = mRatio(1 : size(top, 1) - 2);
    
end