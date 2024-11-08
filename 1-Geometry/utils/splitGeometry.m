function [top, bot] = splitGeometry(xyvals, thres)

    % Split a closed contour into a top half and a bottom half
    [minX, minI] = min(xyvals(:, 1));
    [maxX, maxI] = max(xyvals(:, 1));
    if minI < maxI
        piece1 = xyvals(minI : maxI, :);
        piece2 = [xyvals(maxI : end, :); xyvals(1 : minI, :)];
        piece2 = flip(piece2);
    else
        piece1 = xyvals(maxI : minI, :);
        piece1 = flip(piece1);
        piece2 = [xyvals(minI : end, :); xyvals(1 : maxI, :)];
    end
    while piece1(1, 1) <= minX + thres
        piece1 = piece1(2 : end, :);
    end
    while piece1(end, 1) >= maxX - thres
        piece1 = piece1(1 : end - 1, :);
    end
    while piece2(1, 1) <= minX + thres
        piece2 = piece2(2 : end, :);
    end
    while piece2(end, 1) >= maxX - thres
        piece2 = piece2(1 : end - 1, :);
    end
    if mean(piece1(:, 2)) > mean(piece2(:, 2))
        top = piece1;
        bot = piece2;
    else
        top = piece2;
        bot = piece1;
    end
    while bot(1, 2) >= top(1, 2)
        bot = bot(2 : end, :);
    end
    while bot(end, 2) >= top(end, 2)
        bot = bot(1 : end - 1, :);
    end

end