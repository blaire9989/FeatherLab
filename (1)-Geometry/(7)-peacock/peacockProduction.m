mkdir("data7");

% Color 1
name = "peacock1Cyan";
mkdir("data7/" + name);
mkdir("data7/" + name + "/coefs");
mkdir("data7/" + name + "/geometry");
mkdir("data7/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef1Cyan/coef" + num2str(i) + ".txt");
    id = fopen("data7/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = peacock(40.0, 0.03);
    writematrix(xyvals, "data7/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data7/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data7/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(1, xyvals, info);

% Color 2
name = "peacock2Yellow";
mkdir("data7/" + name);
mkdir("data7/" + name + "/coefs");
mkdir("data7/" + name + "/geometry");
mkdir("data7/" + name + "/render");
for i = 400 : 6 : 694
    coef = readmatrix("coef2Yellow/coef" + num2str(i) + ".txt");
    id = fopen("data7/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = peacock(40.0, 0.03);
    writematrix(xyvals, "data7/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data7/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data7/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(2, xyvals, info);