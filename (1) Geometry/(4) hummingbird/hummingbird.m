function [xyvals, info] = hummingbird(width, cortex, tTop, tReg, tBubble, tSpace, tMembrane, meanMelanin, nL, d)
addpath('../utils/');

    % Compute melanin layer thicknesses and air bubble thicknesses
    tWall = (tReg - tBubble) / 2;
    if tTop == 0
        airThickness = zeros(1, nL);
        for i = 1 : nL
            airThickness(i) = 0.95 * tBubble + 0.10 * tBubble * rand;
        end
        keratinThickness = zeros(1, nL - 1);
        for i = 1 : nL - 1
            keratinThickness(i) = 0.95 * tSpace + 0.10 * tSpace * rand;
        end
        thickness = 2 * cortex + 2 * nL * tWall + sum(airThickness) + sum(keratinThickness);
    else
        b0 = tTop - 2 * tWall;
        if b0 < 0.02
            b0 = b0 + 0.01;
        end
        airThickness = zeros(1, nL);
        airThickness(1) = b0;
        airThickness(2) = 0.20 * b0 + 0.80 * tBubble;
        for i = 3 : nL
            airThickness(i) = 0.95 * tBubble + 0.10 * tBubble * rand;
        end
        keratinThickness = zeros(1, nL - 1);
        keratinThickness(1) = 0.90 * tSpace;
        keratinThickness(2) = 0.95 * tSpace;
        for i = 3 : nL - 1
            keratinThickness(i) = 0.95 * tSpace + 0.10 * tSpace * rand;
        end
        thickness = 2 * cortex + 2 * nL * tWall + sum(airThickness) + sum(keratinThickness);
    end

    % Progressively generate each layer in the barbule
    [yvec1, yvec2, xtop, ytop, xbot, ybot] = generateTopAndBottom(width, thickness, d);
    xyvals = [xtop ytop];
    info = [0 size(xyvals, 1)];
    h1 = cortex;
    for i = 1 : nL
        h2 = h1 + 2 * tWall + airThickness(i);
        [xleft, xright, topLayer, botLayer] = computeLayerExtent(width, thickness, h1, h2, yvec1, yvec2, xtop, ytop, xbot, ybot);
        [xyvals, info] = generateLayer(i, xyvals, info, width, xleft, xright, topLayer, botLayer, tWall, airThickness(i), tMembrane, meanMelanin, d);
        if i < nL
            h1 = h2 + keratinThickness(i);
        end
    end
    info = [info; size(xyvals, 1) size(xbot, 1)];
    xyvals = [xyvals; xbot ybot];

end

function [yvec1, yvec2, xtop, ytop, xbot, ybot] = generateTopAndBottom(width, thickness, d)

    % Generate the barbule top and bottom surfaces, ensuring the structural coloration is accurately preserved
    N = ceil(width / d);
    xvec = transpose(linspace(-0.5 * width, 0.5 * width, N + 1));
    yvec1 = roughness(1.5, -0.5 * width, 0.5 * width, -0.25, 0.25, N + 1, true);
    yvec2 = yvec1 - thickness + roughness(2.0, -0.5 * width, 0.5 * width, -0.05, 0.05, N + 1, true);
    buffer = 1.0;
    d0 = width / N;
    indexL = round(buffer / d0) + 1;
    indexR = round((width - buffer) / d0) + 1;
    xtop = xvec(indexL : indexR);
    ytop = yvec1(indexL : indexR);
    xbot = xtop;
    ybot = yvec2(indexL : indexR);
    angles = linspace(0, 2.0 * pi, 1801);
    slopes = zeros(1, 1800);
    for i = 1 : 1800
        slopes(i) = -cos(angles(i)) / sin(angles(i));
    end

    % Generate the rounded barbule top left corner
    k = (yvec1(indexL + 1) - yvec1(indexL)) / (xvec(indexL + 1) - xvec(indexL));
    [~, index] = min(abs(slopes - k));
    phi0 = angles(index);
    if phi0 > pi
        phi0 = phi0 - pi;
    end
    scale = (xvec(indexL) + 0.5 * width) / (cos(phi0) - cos(2.0 * pi / 3));
    translate = [xvec(indexL) yvec1(indexL)] - scale * [cos(phi0) sin(phi0)];
    n = ceil(scale * abs(phi0 - 2.0 * pi / 3) / d);
    phi = linspace(2.0 * pi / 3, phi0, n + 1);
    for i = n : - 1 : 1
        xtop = [scale * cos(phi(i)) + translate(1); xtop];
        ytop = [scale * sin(phi(i)) + translate(2); ytop];
    end

    % Generate the rounded barbule top right corner
    k = (yvec1(indexR) - yvec1(indexR - 1)) / (xvec(indexR) - xvec(indexR - 1));
    [~, index] = min(abs(slopes - k));
    phi0 = angles(index);
    if phi0 > pi
        phi0 = phi0 - pi;
    end
    scale = (0.5 * width - xvec(indexR)) / (cos(pi / 3) - cos(phi0));
    translate = [xvec(indexR) yvec1(indexR)] - scale * [cos(phi0) sin(phi0)];
    n = ceil(scale * abs(phi0 - pi / 3) / d);
    phi = linspace(phi0, pi / 3, n + 1);
    for i = 2 : n + 1
        xtop = [xtop; scale * cos(phi(i)) + translate(1)];
        ytop = [ytop; scale * sin(phi(i)) + translate(2)];
    end

    % Generate the rounded barbule bottom left corner
    k = (yvec2(indexL + 1) - yvec2(indexL)) / (xvec(indexL + 1) - xvec(indexL));
    [~, index] = min(abs(slopes - k));
    phi0 = angles(index);
    if phi0 < pi
        phi0 = phi0 + pi;
    end
    scale = (xvec(indexL) + 0.5 * width) / (cos(phi0) - cos(4 * pi / 3));
    translate = [xvec(indexL) yvec2(indexL)] - scale * [cos(phi0) sin(phi0)];
    n = ceil(scale * abs(phi0 - 4 * pi / 3) / d);
    phi = linspace(4 * pi / 3, phi0, n + 1);
    for i = n : - 1 : 1
        xbot = [scale * cos(phi(i)) + translate(1); xbot];
        ybot = [scale * sin(phi(i)) + translate(2); ybot];
    end

    % Generate the rounded barbule bottom right corner
    k = (yvec2(indexR) - yvec2(indexR - 1)) / (xvec(indexR) - xvec(indexR - 1));
    [~, index] = min(abs(slopes - k));
    phi0 = angles(index);
    if phi0 < pi
        phi0 = phi0 + pi;
    end
    scale = (0.5 * width - xvec(indexR)) / (cos(5 * pi / 3) - cos(phi0));
    translate = [xvec(indexR) yvec2(indexR)] - scale * [cos(phi0) sin(phi0)];
    n = ceil(scale * abs(phi0 - 5 * pi / 3) / d);
    phi = linspace(phi0, 5 * pi / 3, n + 1);
    for i = 2 : n + 1
        xbot = [xbot; scale * cos(phi(i)) + translate(1)];
        ybot = [ybot; scale * sin(phi(i)) + translate(2)];
    end

    % Generate layers later used for determining layer boundaries in the middle
    for i = 1 : indexL - 1
        count = 1;
        while xtop(count) <= xvec(i)
            count = count + 1;
        end
        if count == 1
            count = 2;
        end
        r = (xvec(i) - xtop(count - 1)) / (xtop(count) - xtop(count - 1));
        yvec1(i) = (1 - r) * ytop(count - 1) + r * ytop(count);
        while xbot(count) <= xvec(i)
            count = count + 1;
        end
        if count == 1
            count = 2;
        end
        r = (xvec(i) - xbot(count - 1)) / (xbot(count) - xbot(count - 1));
        yvec2(i) = (1 - r) * ybot(count - 1) + r * ybot(count);
    end
    for i = indexR + 1 : length(xvec)
        count = 1;
        while count <= length(xtop) && xtop(count) <= xvec(i)
            count = count + 1;
        end
        if count > length(xtop)
            count = length(xtop);
        end
        r = (xvec(i) - xtop(count - 1)) / (xtop(count) - xtop(count - 1));
        yvec1(i) = (1 - r) * ytop(count - 1) + r * ytop(count);
        while count <= length(xbot) && xbot(count) <= xvec(i)
            count = count + 1;
        end
        if count > length(xbot)
            count = length(xbot);
        end
        r = (xvec(i) - xbot(count - 1)) / (xbot(count) - xbot(count - 1));
        yvec2(i) = (1 - r) * ybot(count - 1) + r * ybot(count);
    end

end

function [xleft, xright, topLayer, botLayer] = computeLayerExtent(width, thickness, h1, h2, yvec1, yvec2, xtop, ytop, xbot, ybot)

    % Determine some constrains on the pancake-shaped melanosome placement
    topLayer = (1 - h1 / thickness) * yvec1 + h1 / thickness * yvec2;
    botLayer = (1 - h2 / thickness) * yvec1 + h2 / thickness * yvec2;
    y1 = mean(topLayer);
    y2 = mean(botLayer);

    % Constraints from the layer top boundary
    if y1 <= ytop(1)
        left1 = -0.5 * width;
    else
        count = 1;
        while ytop(count) < y1
            count = count + 1;
        end
        left1 = xtop(count);
    end
    if y1 <= ytop(end)
        right1 = 0.5 * width;
    else
        count = length(ytop);
        while ytop(count) < y1
            count = count - 1;
        end
        right1 = xtop(count);
    end

    % Constraints from the layer bottom boundary
    if y2 >= ybot(1)
        left2 = -0.5 * width;
    else
        count = 1;
        while ybot(count) > y2
            count = count + 1;
        end
        left2 = xbot(count);
    end
    if y2 >= ybot(end)
        right2 = 0.5 * width;
    else
        count = length(ybot);
        while ybot(count) > y2
            count = count - 1;
        end
        right2 = xbot(count);
    end
    xleft = max(left1, left2);
    xright = min(right1, right2);

end

function [xyvals, info] = generateLayer(index, xyvals, info, width, xleft, xright, topLayer, botLayer, tWall, tBubble, tMembrane, meanMelanin, d)

    % Determine the total length of melanosomes
    if index == 1
        ratio = 0.85;
    else
        ratio = 0.90;
    end
    totalMelanin = (xright - xleft) * ratio;

    % Determine the length of each melanosome
    nMelanosome = round(totalMelanin / meanMelanin);
    melaninSize = zeros(nMelanosome, 1);
    for i = 1 : nMelanosome
        melaninSize(i) = 2 * meanMelanin / 3 + 2 * meanMelanin / 3 * rand;
    end
    melaninSize = melaninSize / sum(melaninSize) * totalMelanin;

    % Determine the length of each keratin gap
    gapL = 0.08 + 0.04 * rand;
    gapR = 0.08 + 0.04 * rand;
    totalKeratin = xright - xleft - gapL - gapR - totalMelanin;
    meanKeratin = totalKeratin / (nMelanosome - 1);
    keratinSize = zeros(nMelanosome - 1, 1);
    for i = 1 : nMelanosome - 1
        keratinSize(i) = 2 * meanKeratin / 3 + 2 * meanKeratin / 3 * rand;
    end
    keratinSize = keratinSize / sum(keratinSize) * totalKeratin;

    % Draw the air-filled melanosomes
    N = length(topLayer) - 1;
    d0 = width / N;
    mleft = xleft + gapL;
    for i = 1 : nMelanosome
        % Generate the pancake-shaped melanosome
        mright = mleft + melaninSize(i);
        indexL = round((mleft + 0.5 * width) / d0) + 1;
        indexR = round((mright + 0.5 * width) / d0) + 1;
        h0 = mean(topLayer(indexL : indexR) - botLayer(indexL : indexR));
        indexL = round((mleft + 0.5 * h0 + 0.5 * width) / d0) + 1;
        indexR = round((mright - 0.5 * h0 + 0.5 * width) / d0) + 1;
        melaninBox = makeMelaninContour(topLayer, botLayer, indexL, indexR, width, d);
        
        % Roughly determine the air bubble sizes
        meanAir = tBubble;
        nAir = round((max(melaninBox(:, 1)) - min(melaninBox(:, 1)) - 2 * tWall + tMembrane) / (meanAir + tMembrane));
        meanAir = (max(melaninBox(:, 1)) - min(melaninBox(:, 1)) - 2 * tWall - (nAir - 1) * tMembrane) / nAir;
        airSize = zeros(nAir, 1);
        for j = 1 : nAir
            airSize(j) = 0.7 * meanAir + 0.6 * meanAir * rand;
        end
        airSize = airSize / sum(airSize) * (max(melaninBox(:, 1)) - min(melaninBox(:, 1)) - 2 * tWall - (nAir - 1) * tMembrane);
        if tBubble > 0.03
            bmin = tBubble - 0.01;
            bmax = tBubble + 0.01;
        else
            bmin = tBubble;
            bmax = tBubble;
        end

        % Generate air bubbles
        x1 = min(melaninBox(:, 1)) + tWall;
        xc = zeros(nAir, 1);
        yc = zeros(nAir, 1);
        axis1 = zeros(nAir, 1);
        axis2 = zeros(nAir, 1);
        for j = 1 : nAir
            x2 = x1 + airSize(j);
            x0 = 0.5 * (x1 + x2);
            a1 = 0.5 * airSize(j);
            a2 = 0.5 * (bmin + (bmax - bmin) * rand);
            y0 = [];
            for k = 1 : size(melaninBox, 1)
                next = k + 1;
                if next > size(melaninBox, 1)
                    next = 1;
                end
                if x0 == melaninBox(k, 1)
                    y0 = [y0; melaninBox(k, 2)];
                elseif (melaninBox(k, 1) - x0) * (melaninBox(next, 1) - x0) < 0
                    r = (x0 - melaninBox(k, 1)) / (melaninBox(next, 1) - melaninBox(k, 1));
                    y0 = [y0; (1 - r) * melaninBox(k, 2) + r * melaninBox(next, 2)];
                end
            end
            xc(j) = x0;
            yc(j) = mean(y0);
            axis1(j) = a1;
            axis2(j) = a2;
            x1 = x2 + tMembrane;
        end
        airBubbles = connectAirBubbles(xc, yc, axis1, axis2, pi / 12, d);

        % Put together the layer model
        [melanin1, melanin2] = splitGeometry(melaninBox, 0.01);
        [bubble1, bubble2] = splitGeometry(airBubbles, 0.01);
        info = [info; size(xyvals, 1) size(melanin1, 1)];
        xyvals = [xyvals; melanin1];
        info = [info; size(xyvals, 1) size(bubble1, 1)];
        xyvals = [xyvals; bubble1];
        info = [info; size(xyvals, 1) size(bubble2, 1)];
        xyvals = [xyvals; bubble2];
        info = [info; size(xyvals, 1) size(melanin2, 1)];
        xyvals = [xyvals; melanin2];

        % Update the position pointer
        if  i < nMelanosome
            mleft = mright + keratinSize(i);
        end
    end

end

function melaninBox = makeMelaninContour(topLayer, botLayer, indexL, indexR, width, d)

    % Determine the size of the pancake-shaped melanosome
    N = length(topLayer) - 1;
    xvec = transpose(linspace(-0.5 * width, 0.5 * width, N + 1));
    xtop = xvec(indexL : indexR);
    ytop = topLayer(indexL : indexR);
    xbot = xtop;
    ybot = botLayer(indexL : indexR);
    angles = linspace(0, 2.0 * pi, 1801);
    slopes = zeros(1, 1800);
    for i = 1 : 1800
        slopes(i) = -cos(angles(i)) / sin(angles(i));
    end

    % Model the left end as a circle
    slope1 = (ytop(2) - ytop(1)) / (xtop(2) - xtop(1));
    slope2 = (ybot(2) - ybot(1)) / (xbot(2) - xbot(1));
    [~, index1] = min(abs(slopes - slope1));
    phi1 = angles(index1);
    if phi1 > pi
        phi1 = phi1 - pi;
    end
    pt1 = [cos(phi1) sin(phi1)];
    [~, index2] = min(abs(slopes - slope2));
    phi2 = angles(index2);
    if phi2 < pi
        phi2 = phi2 + pi;
    end
    pt2 = [cos(phi2) sin(phi2)];

    % Scale and translate the left ellipse
    left = [];
    if pt1(1) <= pt2(1)
        % Add flat segment at the bottom
        c1 = (ybot(2) - ybot(1)) * (pt1(1) - pt2(1)) - (xbot(2) - xbot(1)) * (pt1(2) - pt2(2));
        c2 = (ybot(2) - ybot(1)) * (xbot(1) - xtop(1)) - (xbot(2) - xbot(1)) * (ybot(1) - ytop(1));
        scale = -c2 / c1;
        translate = [xtop(1) ytop(1)] - scale * pt1;
        n = ceil(scale * abs(phi1 - phi2) / d);
        phi = linspace(phi1, phi2, n + 1);
        for i = 2 : n
            left = [left; scale * [cos(phi(i)) sin(phi(i))] + translate];
        end
    else
        % Add flat segment at the top
        c1 = (ytop(2) - ytop(1)) * (pt1(1) - pt2(1)) - (xtop(2) - xtop(1)) * (pt1(2) - pt2(2));
        c2 = (ytop(2) - ytop(1)) * (xtop(1) - xbot(1)) - (xtop(2) - xtop(1)) * (ytop(1) - ybot(1));
        scale = c2 / c1;
        translate = [xbot(1) ybot(1)] - scale * pt2;
        n = ceil(scale * abs(phi1 - phi2) / d);
        phi = linspace(phi1, phi2, n + 1);
        for i = 2 : n
            left = [left; scale * [cos(phi(i)) sin(phi(i))] + translate];
        end
    end

    % Model the right end as a circle
    slope1 = (ytop(end) - ytop(end - 1)) / (xtop(end) - xtop(end - 1));
    slope2 = (ybot(end) - ybot(end - 1)) / (xbot(end) - xbot(end - 1));
    [~, index1] = min(abs(slopes - slope1));
    phi1 = angles(index1);
    if phi1 > pi
        phi1 = phi1 - pi;
    end
    pt1 = [cos(phi1) sin(phi1)];
    phi1 = phi1 + 2.0 * pi;
    [~, index2] = min(abs(slopes - slope2));
    phi2 = angles(index2);
    if phi2 < pi
        phi2 = phi2 + pi;
    end
    pt2 = [cos(phi2) sin(phi2)];

    % Scale and translate the right ellipse
    right = [];
    if pt2(1) <= pt1(1)
        % Add flat segment at the bottom
        c1 = (ybot(end - 1) - ybot(end)) * (pt1(1) - pt2(1)) - (xbot(end - 1) - xbot(end)) * (pt1(2) - pt2(2));
        c2 = (ybot(end - 1) - ybot(end)) * (xbot(end) - xtop(end)) - (xbot(end - 1) - xbot(end)) * (ybot(end) - ytop(end));
        scale = -c2 / c1;
        translate = [xtop(end) ytop(end)] - scale * pt1;
        n = ceil(scale * abs(phi1 - phi2) / d);
        phi = linspace(phi2, phi1, n + 1);
        for i = 2 : n
            right = [right; scale * [cos(phi(i)) sin(phi(i))] + translate];
        end
    else
        % Add flat segment at the top
        c1 = (ytop(end - 1) - ytop(end)) * (pt1(1) - pt2(1)) - (xtop(end - 1) - xtop(end)) * (pt1(2) - pt2(2));
        c2 = (ytop(end - 1) - ytop(end)) * (xtop(end) - xbot(end)) - (xtop(end - 1) - xtop(end)) * (ytop(end) - ybot(end));
        scale = c2 / c1;
        translate = [xbot(end) ybot(end)] - scale * pt2;
        n = ceil(scale * abs(phi1 - phi2) / d);
        phi = linspace(phi2, phi1, n + 1);
        for i = 2 : n
            right = [right; scale * [cos(phi(i)) sin(phi(i))] + translate];
        end
    end
    melaninBox = [flip(xtop) flip(ytop); left; xbot ybot; right];

end

function airBubbles = connectAirBubbles(xc, yc, axis1, axis2, critical, d)

    % Put in the leftmost bubble
    airBubbles = [];
    phi1 = 2.0 * pi - critical;
    phi2 = critical;
    n = max(7, ceil(max(axis1(1), axis2(1)) * (phi1 - phi2) / d));
    phi = linspace(phi1, phi2, n + 1);
    for i = 1 : n
        airBubbles = [airBubbles; xc(1) + axis1(1) * cos(phi(i)) yc(1) + axis2(1) * sin(phi(i))];
    end

    % Put in the top boundary
    for i = 1 : length(xc) - 1
        pt1 = [xc(i) + axis1(i) * cos(critical) yc(i) + axis2(i) * sin(critical)];
        pt2 = [xc(i + 1) + axis1(i + 1) * cos(pi - critical) yc(i + 1) + axis2(i + 1) * sin(pi - critical)];
        dist = sqrt((pt1(1) - pt2(1))^2 + (pt1(2) - pt2(2))^2);
        n = ceil(dist / d);
        xvec = linspace(pt1(1), pt2(1), n + 1);
        yvec = linspace(pt1(2), pt2(2), n + 1);
        airBubbles = [airBubbles; transpose(xvec(1 : end - 1)) transpose(yvec(1 : end - 1))];
        if i < length(xc) - 1
            phi1 = pi - critical;
            phi2 = critical;
            n = max(3, ceil(max(axis1(i + 1), axis2(i + 1)) * (phi1 - phi2) / d));
            phi = linspace(phi1, phi2, n + 1);
            for j = 1 : n
                airBubbles = [airBubbles; xc(i + 1) + axis1(i + 1) * cos(phi(j)) yc(i + 1) + axis2(i + 1) * sin(phi(j))];
            end
        end
    end
    
    % Put in the rightmost bubble
    phi1 = pi - critical;
    phi2 = critical - pi;
    n = max(7, ceil(max(axis1(end), axis2(end)) * (phi1 - phi2) / d));
    phi = linspace(phi1, phi2, n + 1);
    for i = 1 : n
        airBubbles = [airBubbles; xc(end) + axis1(end) * cos(phi(i)) yc(end) + axis2(end) * sin(phi(i))];
    end

    % Put in the bottom boundary
    for i = length(xc) : -1 : 2
        pt1 = [xc(i) + axis1(i) * cos(critical - pi) yc(i) + axis2(i) * sin(critical - pi)];
        pt2 = [xc(i - 1) + axis1(i - 1) * cos(-critical) yc(i - 1) + axis2(i - 1) * sin(-critical)];
        dist = sqrt((pt1(1) - pt2(1))^2 + (pt1(2) - pt2(2))^2);
        n = ceil(dist / d);
        xvec = linspace(pt1(1), pt2(1), n + 1);
        yvec = linspace(pt1(2), pt2(2), n + 1);
        airBubbles = [airBubbles; transpose(xvec(1 : end - 1)) transpose(yvec(1 : end - 1))];
        if i > 2
            phi1 = -critical;
            phi2 = critical - pi;
            n = max(3, ceil(max(axis1(i - 1), axis2(i - 1)) * (phi1 - phi2) / d));
            phi = linspace(phi1, phi2, n + 1);
            for j = 1 : n
                airBubbles = [airBubbles; xc(i - 1) + axis1(i - 1) * cos(phi(j)) yc(i - 1) + axis2(i - 1) * sin(phi(j))];
            end
        end
    end

end