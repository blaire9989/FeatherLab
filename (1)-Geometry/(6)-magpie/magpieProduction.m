mkdir("data6");

% Color 1
name = "magpie1Forest";
mkdir("data6/" + name);
mkdir("data6/" + name + "/coefs");
mkdir("data6/" + name + "/geometry");
mkdir("data6/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef1Forest/coef" + num2str(i) + ".txt");
    id = fopen("data6/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = magpie(40.0, 3.6, 1.0, 0.03);
    writematrix(xyvals, "data6/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data6/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data6/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(1, xyvals, info);

% Color 2
name = "magpie2Lime";
mkdir("data6/" + name);
mkdir("data6/" + name + "/coefs");
mkdir("data6/" + name + "/geometry");
mkdir("data6/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef2Lime/coef" + num2str(i) + ".txt");
    id = fopen("data6/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = magpie(40.0, 3.6, 1.0, 0.03);
    writematrix(xyvals, "data6/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data6/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data6/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(2, xyvals, info);

% Color 3
name = "magpie3Yellow";
mkdir("data6/" + name);
mkdir("data6/" + name + "/coefs");
mkdir("data6/" + name + "/geometry");
mkdir("data6/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef3Yellow/coef" + num2str(i) + ".txt");
    id = fopen("data6/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = magpie(40.0, 3.6, 1.0, 0.03);
    writematrix(xyvals, "data6/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data6/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data6/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(3, xyvals, info);

% Color 4
name = "magpie4Rose";
mkdir("data6/" + name);
mkdir("data6/" + name + "/coefs");
mkdir("data6/" + name + "/geometry");
mkdir("data6/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef4Rose/coef" + num2str(i) + ".txt");
    id = fopen("data6/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = magpie(40.0, 3.3, 0.9, 0.03);
    writematrix(xyvals, "data6/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data6/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data6/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(4, xyvals, info);

% Color 5
name = "magpie5Purple";
mkdir("data6/" + name);
mkdir("data6/" + name + "/coefs");
mkdir("data6/" + name + "/geometry");
mkdir("data6/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef5Purple/coef" + num2str(i) + ".txt");
    id = fopen("data6/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = magpie(40.0, 3.3, 0.8, 0.03);
    writematrix(xyvals, "data6/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data6/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data6/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(5, xyvals, info);

% Color 6
name = "magpie6Violet";
mkdir("data6/" + name);
mkdir("data6/" + name + "/coefs");
mkdir("data6/" + name + "/geometry");
mkdir("data6/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef6Violet/coef" + num2str(i) + ".txt");
    id = fopen("data6/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = magpie(40.0, 3.3, 0.7, 0.03);
    writematrix(xyvals, "data6/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data6/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data6/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(6, xyvals, info);

% Color 7
name = "magpie7Blue";
mkdir("data6/" + name);
mkdir("data6/" + name + "/coefs");
mkdir("data6/" + name + "/geometry");
mkdir("data6/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef7Blue/coef" + num2str(i) + ".txt");
    id = fopen("data6/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = magpie(40.0, 3.0, 0.6, 0.03);
    writematrix(xyvals, "data6/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data6/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data6/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(7, xyvals, info);

% Color 8
name = "magpie8Teal";
mkdir("data6/" + name);
mkdir("data6/" + name + "/coefs");
mkdir("data6/" + name + "/geometry");
mkdir("data6/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef8Teal/coef" + num2str(i) + ".txt");
    id = fopen("data6/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = magpie(40.0, 3.0, 0.5, 0.03);
    writematrix(xyvals, "data6/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data6/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data6/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(8, xyvals, info);