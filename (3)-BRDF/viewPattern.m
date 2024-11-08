function colors = viewPattern(type, name, index, D65)

    % Read the output data files describing scattering at 50 wavelengths
    resolution = 400;
    spec1 = zeros(resolution, 50);
    spec2 = zeros(resolution, 50);
    spec3 = zeros(resolution, 50);
    spec4 = zeros(resolution, 50);
    spec5 = zeros(resolution, 50);
    for i = 1 : 50
        wvl = 400 + 6 * (i - 1);
        id = fopen("data" + num2str(type) + "/" + name + "/visual/pattern_" + num2str(wvl) + "_" + num2str(index) + ".binary", 'r');
        pattern = fread(id, [resolution 5], 'float');
        fclose(id);
        spec1(:, i) = pattern(:, 1);
        spec2(:, i) = pattern(:, 2);
        spec3(:, i) = pattern(:, 3);
        spec4(:, i) = pattern(:, 4);
        spec5(:, i) = pattern(:, 5);
    end

    % Convert the spectral data into RGB
    spec1 = ToXYZ(spec1, 50, 400, 700, D65);
    spec1 = XYZToRGB(spec1);
    spec2 = ToXYZ(spec2, 50, 400, 700, D65);
    spec2 = XYZToRGB(spec2);
    spec3 = ToXYZ(spec3, 50, 400, 700, D65);
    spec3 = XYZToRGB(spec3);
    spec4 = ToXYZ(spec4, 50, 400, 700, D65);
    spec4 = XYZToRGB(spec4);
    spec5 = ToXYZ(spec5, 50, 400, 700, D65);
    spec5 = XYZToRGB(spec5);
    slice1 = zeros(resolution / 5, resolution, 3);
    slice2 = zeros(resolution / 5, resolution, 3);
    slice3 = zeros(resolution / 5, resolution, 3);
    slice4 = zeros(resolution / 5, resolution, 3);
    slice5 = zeros(resolution / 5, resolution, 3);
    for i = 1 : resolution / 5
        for j = 1 : resolution
            for k = 1 : 3
                slice1(i, j, k) = spec1(j, k);
                slice2(i, j, k) = spec2(j, k);
                slice3(i, j, k) = spec3(j, k);
                slice4(i, j, k) = spec4(j, k);
                slice5(i, j, k) = spec5(j, k);
            end
        end
    end
    white1 = ones(10, resolution, 3);
    colors = [white1; slice1; white1; slice2; white1; slice3; white1; slice4; white1; slice5; white1];
    white2 = ones(size(colors, 1), 10, 3);
    colors = [white2 colors white2];
    colors(colors < 0) = 0;
    colors = colors .^ (1 / 2.2);

end