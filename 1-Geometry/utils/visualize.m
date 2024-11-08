function visualize(fig_num, xyvals, info)
    
    % Visualize a multi-component barbule structure
    figure(fig_num);
    for i = 1 : size(info, 1)
        start = info(i, 1);
        total = info(i, 2);
        plot(xyvals(start + 1 : start + total, 1), xyvals(start + 1 : start + total, 2)); hold on; axis equal;
    end

end