function [xyvals, info, mRatio] = bronzewing(width, cortex, Dm, spacing, nL, d)
addpath('../utils/');

    % Generate the outer contour
    N = round(width / d);
    if rem(N, 2) == 1
        N = N + 1;
    end
    period = Dm + spacing;
    thickness = 2 * cortex + 2 * nL * period + 2.5 * nL * period;
    xtop = transpose(linspace(-0.5 * width, 0.5 * width, N + 1));
    ytop = zeros(N + 1, 1);
    xbot = transpose(linspace(-0.5 * width, 0.5 * width, N + 1));
    ybot = -thickness * ones(N + 1, 1);
    contour = fullContour(xtop, ytop, xbot, ybot, d, 2.0 + 1.0 * rand);
    perturb = roughness(1.5, -width, width, -0.2, 0.2, size(contour, 1), false);
    contour(:, 2) = contour(:, 2) + perturb(1 : size(contour, 1));
    [top, bot] = splitGeometry(contour, 1.0);
    bot = bot * width / (max(top(:, 1)) - min(top(:, 1)));
    top = top * width / (max(top(:, 1)) - min(top(:, 1)));
    xshift = 0.5 * (min(top(:, 1)) + max(top(:, 1)));
    yshift = mean(top(:, 2));
    top(:, 1) = top(:, 1) - xshift;
    top(:, 2) = top(:, 2) - yshift;
    bot(:, 1) = bot(:, 1) - xshift;
    bot(:, 2) = bot(:, 2) - mean(bot(:, 2)) - thickness;
    layers = cell(4 * nL, 1);

    % Arrange melanosome layers close to the barbule top surface
    guide1 = top;
    slopeL = (guide1(2, 2) - guide1(1, 2)) / (guide1(2, 1) - guide1(1, 1));
    x0 = guide1(1, 1) - nL * 1.0;
    y0 = guide1(1, 2) - nL * 1.0 * slopeL;
    n = ceil(nL * 1.0 * sqrt(1 + slopeL * slopeL) / d);
    xExtra = linspace(x0, guide1(1, 1), n + 1);
    yExtra = linspace(y0, guide1(1, 2), n + 1);
    guide1 = [transpose(xExtra(1 : n)) transpose(yExtra(1 : n)); guide1];
    slopeR = (guide1(end, 2) - guide1(end - 1, 2)) / (guide1(end, 1) - guide1(end - 1, 1));
    x0 = guide1(end, 1) + nL * 1.0;
    y0 = guide1(end, 2) + nL * 1.0 * slopeR;
    n = ceil(nL * 1.0 * sqrt(1 + slopeR * slopeR) / d);
    xExtra = linspace(guide1(end, 1), x0, n + 1);
    yExtra = linspace(guide1(end, 2), y0, n + 1);
    guide1 = [guide1; transpose(xExtra(2 : n + 1)) transpose(yExtra(2 : n + 1))];
    for i = 1 : nL
        guide = zeros(size(guide1, 1), 2);
        for j = 1 : size(guide1, 1)
            if j == 1
                tvec = [guide1(2, 1) - guide1(1, 1) guide1(2, 2) - guide1(1, 2)];
            elseif j == size(guide1, 1)
                tvec = [guide1(end, 1) - guide1(end - 1, 1) guide1(end, 2) - guide1(end - 1, 2)];
            else
                tvec = [guide1(j + 1, 1) - guide1(j - 1, 1) guide1(j + 1, 2) - guide1(j - 1, 2)];
            end
            tvec = tvec / norm(tvec);
            guide(j, :) = guide1(j, :) - (cortex + 0.5 * Dm + (i - 1) * period) * [-tvec(2) tvec(1)];
        end
        [xc, yc, vD] = growInterior(guide, Dm, 0.1);
        melanin = generateLayer(xc, yc, vD, d);
        [mel1, mel2] = splitGeometry(melanin, 0.1);
        index1 = 1;
        while mel1(index1, 1) < -0.5 * width
            index1 = index1 + 1;
        end
        index2 = size(mel1, 1);
        while mel1(index2, 1) > 0.5 * width
            index2 = index2 - 1;
        end
        mel1 = mel1(index1 : index2, :);
        index1 = 1;
        while mel2(index1, 1) < -0.5 * width
            index1 = index1 + 1;
        end
        index2 = size(mel2, 1);
        while mel2(index2, 1) > 0.5 * width
            index2 = index2 - 1;
        end
        mel2 = mel2(index1 : index2, :);
        layers{2 * i - 1} = mel1;
        layers{2 * i} = mel2;
    end

    % Arrange melanosome layers close to the barbule bottom surface
    guide2 = bot;
    slopeL = (guide2(2, 2) - guide2(1, 2)) / (guide2(2, 1) - guide2(1, 1));
    x0 = guide2(1, 1) - nL * 1.0;
    y0 = guide2(1, 2) - nL * 1.0 * slopeL;
    n = ceil(nL * 1.0 * sqrt(1 + slopeL * slopeL) / d);
    xExtra = linspace(x0, guide2(1, 1), n + 1);
    yExtra = linspace(y0, guide2(1, 2), n + 1);
    guide2 = [transpose(xExtra(1 : n)) transpose(yExtra(1 : n)); guide2];
    slopeR = (guide2(end, 2) - guide2(end - 1, 2)) / (guide2(end, 1) - guide2(end - 1, 1));
    x0 = guide2(end, 1) + nL * 1.0;
    y0 = guide2(end, 2) + nL * 1.0 * slopeR;
    n = ceil(nL * 1.0 * sqrt(1 + slopeR * slopeR) / d);
    xExtra = linspace(guide2(end, 1), x0, n + 1);
    yExtra = linspace(guide2(end, 2), y0, n + 1);
    guide2 = [guide2; transpose(xExtra(2 : n + 1)) transpose(yExtra(2 : n + 1))];
    for i = 1 : nL
        guide = zeros(size(guide2, 1), 2);
        for j = 1 : size(guide2, 1)
            if j == 1
                tvec = [guide2(2, 1) - guide2(1, 1) guide2(2, 2) - guide2(1, 2)];
            elseif j == size(guide2, 1)
                tvec = [guide2(end, 1) - guide2(end - 1, 1) guide2(end, 2) - guide2(end - 1, 2)];
            else
                tvec = [guide2(j + 1, 1) - guide2(j - 1, 1) guide2(j + 1, 2) - guide2(j - 1, 2)];
            end
            tvec = tvec / norm(tvec);
            guide(j, :) = guide2(j, :) + (cortex + 0.5 * Dm + (nL - i) * period) * [-tvec(2) tvec(1)];
        end
        [xc, yc, vD] = growInterior(guide, Dm, 0.1);
        melanin = generateLayer(xc, yc, vD, d);
        [mel1, mel2] = splitGeometry(melanin, 0.1);
        index1 = 1;
        while mel1(index1, 1) < -0.5 * width
            index1 = index1 + 1;
        end
        index2 = size(mel1, 1);
        while mel1(index2, 1) > 0.5 * width
            index2 = index2 - 1;
        end
        mel1 = mel1(index1 : index2, :);
        index1 = 1;
        while mel2(index1, 1) < -0.5 * width
            index1 = index1 + 1;
        end
        index2 = size(mel2, 1);
        while mel2(index2, 1) > 0.5 * width
            index2 = index2 - 1;
        end
        mel2 = mel2(index1 : index2, :);
        layers{2 * nL + 2 * i - 1} = mel1;
        layers{2 * nL + 2 * i} = mel2;
    end

    % Put together the geometric model
    info = [0 size(top, 1)];
    xyvals = top;
    for i = 1 : 4 * nL
        info = [info; size(xyvals, 1) size(layers{i}, 1)];
        xyvals = [xyvals; layers{i}];
    end
    info = [info; size(xyvals, 1) size(bot, 1)];
    xyvals = [xyvals; bot];

    % Model the percentage of melanin in the interior
    total = size(top, 1) - 2;
    if rem(total, 2) == 1
        total = total + 1;
    end
    mRatio = roughness(1.0, -0.5 * width, 0.5 * width, 0.1, 0.2, total, false);
    mRatio = mRatio(1 : size(top, 1) - 2);

end