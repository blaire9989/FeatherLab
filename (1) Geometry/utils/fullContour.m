function contour = fullContour(xtop, ytop, xbot, ybot, d, longaxis)

    % Connect the top and bottom surfaces into a full barbule contour
    a = longaxis;
    b = 1.0;
    phi = linspace(0, 2 * pi, 1801);
    slopes = zeros(1, 1800);
    for i = 1 : 1800
        sx = -a * sin(phi(i));
        sy = b * cos(phi(i));
        slopes(i) = sy / sx;
    end

    % Model the barbule's left end as an ellipse
    slope1 = (ytop(2) - ytop(1)) / (xtop(2) - xtop(1));
    slope2 = (ybot(2) - ybot(1)) / (xbot(2) - xbot(1));
    [~, index1] = min(abs(slopes - slope1));
    phi1 = phi(index1);
    if phi1 > pi
        phi1 = phi1 - pi;
    end
    pt1 = [a * cos(phi1) b * sin(phi1)];
    [~, index2] = min(abs(slopes - slope2));
    phi2 = phi(index2);
    if phi2 < pi
        phi2 = phi2 + pi;
    end
    pt2 = [a * cos(phi2) b * sin(phi2)];
    L = ellipseArc(a, b, phi1, phi2);

    % Scale and translate the ellipse on the left
    s1 = (pt2(2) - pt1(2)) / (pt2(1) - pt1(1));
    s2 = (ybot(1) - ytop(1)) / (xbot(1) - xtop(1));
    left = [];
    if (s1 > 0 && s2 < 0) || (s1 * s2 > 0 && s1 < s2)
        % Need to add flat segment at the bottom
        c1 = (ybot(2) - ybot(1)) * (pt1(1) - pt2(1)) - (xbot(2) - xbot(1)) * (pt1(2) - pt2(2));
        c2 = (ybot(2) - ybot(1)) * (xbot(1) - xtop(1)) - (xbot(2) - xbot(1)) * (ybot(1) - ytop(1));
        scale = -c2 / c1;
        translate = [xtop(1) ytop(1)] - scale * pt1;
        num = ceil(scale * L / d);
        d0 = scale * L / num;
        count = 1;
        phi0 = phi1;
        while count <= num
            pt = scale * [a * cos(phi0) b * sin(phi0)] + translate;
            left = [left; pt];
            delta_phi = d0 / (scale * sqrt(a^2 * sin(phi0)^2 + b^2 * cos(phi0)^2));
            phi0 = phi0 + delta_phi;
            count = count + 1;
        end
        l1 = scale * [a * cos(phi2) b * sin(phi2)] + translate;
        l2 = [xbot(1) ybot(1)];
        dist = sqrt((l1(1) - l2(1))^2 + (l1(2) - l2(2))^2);
        num = ceil(dist / d);
        xvec = linspace(l1(1), l2(1), num + 1);
        yvec = linspace(l1(2), l2(2), num + 1);
        left = [left(2 : end, :); transpose(xvec(1 : end - 1)) transpose(yvec(1 : end - 1))];
    else
        % Need to add flat segment at the top
        c1 = (ytop(2) - ytop(1)) * (pt1(1) - pt2(1)) - (xtop(2) - xtop(1)) * (pt1(2) - pt2(2));
        c2 = (ytop(2) - ytop(1)) * (xtop(1) - xbot(1)) - (xtop(2) - xtop(1)) * (ytop(1) - ybot(1));
        scale = c2 / c1;
        translate = [xbot(1) ybot(1)] - scale * pt2;
        num = ceil(scale * L / d);
        d0 = scale * L / num;
        count = 1;
        phi0 = phi1;
        while count <= num
            pt = scale * [a * cos(phi0) b * sin(phi0)] + translate;
            left = [left; pt];
            delta_phi = d0 / (scale * sqrt(a^2 * sin(phi0)^2 + b^2 * cos(phi0)^2));
            phi0 = phi0 + delta_phi;
            count = count + 1;
        end
        l1 = [xtop(1) ytop(1)];
        l2 = left(1, :);
        dist = sqrt((l1(1) - l2(1))^2 + (l1(2) - l2(2))^2);
        num = ceil(dist / d);
        xvec = linspace(l1(1), l2(1), num + 1);
        yvec = linspace(l1(2), l2(2), num + 1);
        left = [transpose(xvec(2 : end - 1)) transpose(yvec(2 : end - 1)); left];
    end

    % Model the barbule's right end as an ellipse
    slope1 = (ytop(end) - ytop(end - 1)) / (xtop(end) - xtop(end - 1));
    slope2 = (ybot(end) - ybot(end - 1)) / (xbot(end) - xbot(end - 1));
    [~, index1] = min(abs(slopes - slope1));
    phi1 = phi(index1);
    if phi1 > pi
        phi1 = phi1 - pi;
    end
    pt1 = [a * cos(phi1) b * sin(phi1)];
    phi1 = phi1 + 2 * pi;
    [~, index2] = min(abs(slopes - slope2));
    phi2 = phi(index2);
    if phi2 < pi
        phi2 = phi2 + pi;
    end
    pt2 = [a * cos(phi2) b * sin(phi2)];
    L = ellipseArc(a, b, phi2, phi1);

    % Scale and translate the ellipse on the right
    s1 = (pt2(2) - pt1(2)) / (pt2(1) - pt1(1));
    s2 = (ybot(end) - ytop(end)) / (xbot(end) - xtop(end));
    right = [];
    if (s1 < 0 && s2 > 0) || (s1 * s2 > 0 && s1 > s2)
        % Need to add flat segment at the bottom
        c1 = (ybot(end - 1) - ybot(end)) * (pt1(1) - pt2(1)) - (xbot(end - 1) - xbot(end)) * (pt1(2) - pt2(2));
        c2 = (ybot(end - 1) - ybot(end)) * (xbot(end) - xtop(end)) - (xbot(end - 1) - xbot(end)) * (ybot(end) - ytop(end));
        scale = -c2 / c1;
        translate = [xtop(end) ytop(end)] - scale * pt1;
        num = ceil(scale * L / d);
        d0 = scale * L / num;
        count = 1;
        phi0 = phi2;
        while count <= num
            pt = scale * [a * cos(phi0) b * sin(phi0)] + translate;
            right = [right; pt];
            delta_phi = d0 / (scale * sqrt(a^2 * sin(phi0)^2 + b^2 * cos(phi0)^2));
            phi0 = phi0 + delta_phi;
            count = count + 1;
        end
        l1 = [xbot(end) ybot(end)];
        l2 = right(1, :);
        dist = sqrt((l1(1) - l2(1))^2 + (l1(2) - l2(2))^2);
        num = ceil(dist / d);
        xvec = linspace(l1(1), l2(1), num + 1);
        yvec = linspace(l1(2), l2(2), num + 1);
        right = [transpose(xvec(2 : end - 1)) transpose(yvec(2 : end - 1)); right];
    else
        % Need to add flat segment at the top
        c1 = (ytop(end - 1) - ytop(end)) * (pt1(1) - pt2(1)) - (xtop(end - 1) - xtop(end)) * (pt1(2) - pt2(2));
        c2 = (ytop(end - 1) - ytop(end)) * (xtop(end) - xbot(end)) - (xtop(end - 1) - xtop(end)) * (ytop(end) - ybot(end));
        scale = c2 / c1;
        translate = [xbot(end) ybot(end)] - scale * pt2;
        num = ceil(scale * L / d);
        d0 = scale * L / num;
        count = 1;
        phi0 = phi2;
        while count <= num
            pt = scale * [a * cos(phi0) b * sin(phi0)] + translate;
            right = [right; pt];
            delta_phi = d0 / (scale * sqrt(a^2 * sin(phi0)^2 + b^2 * cos(phi0)^2));
            phi0 = phi0 + delta_phi;
            count = count + 1;
        end
        l1 = scale * [a * cos(phi1) b * sin(phi1)] + translate;
        l2 = [xtop(end) ytop(end)];
        dist = sqrt((l1(1) - l2(1))^2 + (l1(2) - l2(2))^2);
        num = ceil(dist / d);
        xvec = linspace(l1(1), l2(1), num + 1);
        yvec = linspace(l1(2), l2(2), num + 1);
        right = [right(2 : end, :); transpose(xvec(1 : end - 1)) transpose(yvec(1 : end - 1))];
    end
    contour = [flip(xtop) flip(ytop); left; xbot ybot; right];

end