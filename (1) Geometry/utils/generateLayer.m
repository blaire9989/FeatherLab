function layer = generateLayer(xc, yc, vD, d)

    % Generate one layer of interior melanosomes
    anchorTop = zeros(length(vD) - 1, 2);
    anchorBot = zeros(length(vD) - 1, 2);
    for i = 1 : length(vD) - 1
        vec = [xc(i + 1) - xc(i) yc(i + 1) - yc(i)];
        vec = vec / norm(vec);
        pt1 = [xc(i) yc(i)] + 0.5 * vD(i) * vec;
        pt2 = [xc(i + 1) yc(i + 1)] - vD(i + 1) / 2 * vec;
        pt = (pt1 + pt2) / 2;
        perp = [-vec(2) vec(1)];
        anchorTop(i, :) = pt + min(vD(i), vD(i + 1)) / 3 * perp;
        anchorBot(i, :) = pt - min(vD(i), vD(i + 1)) / 3 * perp;
    end
    
    % Build the approximated melanin layer structure
    phi1 = findPhi(xc(end), yc(end), 0.5 * vD(end), anchorTop(end, :));
    phi2 = findPhi(xc(end), yc(end), 0.5 * vD(end), anchorBot(end, :));
    layer = drawArc(xc(end), yc(end), 0.5 * vD(end), d, phi1, phi2);
    for i = length(vD) - 1 : -1 : 2
        phi1 = findPhi(xc(i), yc(i), 0.5 * vD(i), anchorTop(i - 1, :));
        phi2 = findPhi(xc(i), yc(i), 0.5 * vD(i), anchorTop(i, :));
        arc1 = drawArc(xc(i), yc(i), 0.5 * vD(i), d, phi1, phi2);
        phi3 = findPhi(xc(i), yc(i), 0.5 * vD(i), anchorBot(i, :));
        phi4 = findPhi(xc(i), yc(i), 0.5 * vD(i), anchorBot(i - 1, :));
        arc2 = drawArc(xc(i), yc(i), 0.5 * vD(i), d, phi3, phi4);
        layer = [arc1; layer; arc2];
    end
    phi1 = findPhi(xc(1), yc(1), 0.5 * vD(1), anchorBot(1, :));
    phi2 = findPhi(xc(1), yc(1), 0.5 * vD(1), anchorTop(1, :));
    arc = drawArc(xc(1), yc(1), 0.5 * vD(1), d, phi1, phi2);
    layer = [layer; arc];

end

function phi = findPhi(x0, y0, r0, pt)

    % Find the point on the circle that is closest to the given point
    values = linspace(0, 359, 360) * pi / 180;
    xvec = x0 + r0 * cos(values);
    yvec = y0 + r0 * sin(values);
    dist = zeros(360, 1);
    for i = 1 : 360
        dist(i) = sqrt((xvec(i) - pt(1))^2 + (yvec(i) - pt(2))^2);
    end
    [~, index] = min(dist);
    phi = values(index);

end

function arc = drawArc(x0, y0, r0, d, phi1, phi2)

    % Generate a section (arc) from a circle
    if phi1 - phi2 > pi
        phi1 = phi1 - 2 * pi;
    elseif phi1 - phi2 < -pi
        phi1 = phi1 + 2 * pi;
    end
    n = max(4, ceil(abs(phi1 - phi2) * r0 / d));
    phi = linspace(phi1, phi2, n + 1);
    arc = zeros(n, 2);
    for i = 1 : n
        arc(i, 1) = x0 + r0 * cos(phi(i));
        arc(i, 2) = y0 + r0 * sin(phi(i));
    end

end