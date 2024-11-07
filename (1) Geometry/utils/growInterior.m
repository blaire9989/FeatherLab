function [xc, yc, vD] = growInterior(guide, Dm, randomness)

    % Put melanosomes along the guide geometry
    xc = guide(1, 1);
    yc = guide(1, 2);
    vD = (1 - randomness) * Dm + 2 * randomness * Dm * rand;
    curr = 1;
    next = 2;
    while true
        diameter = (1 - randomness) * Dm + 2 * randomness * Dm * rand;
        d1 = vD(end) / 2 + diameter / 2;
        d2 = sqrt((xc(end) - guide(next, 1))^2 + (yc(end) - guide(next, 2))^2);
        if d1 <= d2
            tvec = [guide(next, 1) - guide(curr, 1); guide(next, 2) - guide(curr, 2)];
            tvec = tvec / norm(tvec);
            newXY = [xc(end); yc(end)] + d1 * tvec;
        else
            if next == size(guide, 1)
                break;
            end
            while next <= size(guide, 1) && sqrt((xc(end) - guide(next, 1))^2 + (yc(end) - guide(next, 2))^2) < d1
                next = next + 1;
            end
            if next > size(guide, 1)
                break;
            else
                curr = next - 1;
            end
            d3 = d1 - sqrt((xc(end) - guide(curr, 1))^2 + (yc(end) - guide(curr, 2))^2);
            d4 = sqrt((guide(next, 1) - guide(curr, 1))^2 + (guide(next, 2) - guide(curr, 2))^2);
            if d3 <= d4
                tvec = [guide(next, 1) - guide(curr, 1); guide(next, 2) - guide(curr, 2)];
                tvec = tvec / norm(tvec);
                newXY = [guide(curr, 1); guide(curr, 2)] + d3 * tvec;
            else
                newXY = [guide(next, 1); guide(next, 2)];
            end
        end
        xc = [xc; newXY(1)];
        yc = [yc; newXY(2)];
        vD = [vD; diameter];
    end

end