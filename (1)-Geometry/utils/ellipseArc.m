function L = ellipseArc(a, b, t1, t2)

    % Compute the arc-length of a section from an ellipse
    nseg = 100000;
    nseg = nseg * (t2 - t1) / (2 * pi);
    tt = t1 : (t2 - t1) / nseg : t2;
    xx = a * cos(tt);
    yy = b * sin(tt);
    dx = diff(xx);
    dy = diff(yy);
    L = sum(sqrt(dx.^2 + dy.^2));

end