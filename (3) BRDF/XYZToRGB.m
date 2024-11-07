function rgb = XYZToRGB(xyz)
    rgb = zeros(size(xyz));
    rgb(:,1) = 3.240479 * xyz(:,1) - 1.537150 * xyz(:,2) - 0.498535 * xyz(:, 3);
    rgb(:,2) = -0.969256 * xyz(:,1) + 1.875991 * xyz(:,2) + 0.041556 * xyz(:, 3);
    rgb(:,3) = 0.055648 * xyz(:,1) - 0.204043 * xyz(:,2) + 1.057311 * xyz(:, 3);
end