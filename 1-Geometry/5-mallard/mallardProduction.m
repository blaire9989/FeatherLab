mkdir("data5");

% Color 1
name = "mallard1Lime";
mkdir("data5/" + name);
mkdir("data5/" + name + "/coefs");
mkdir("data5/" + name + "/geometry");
mkdir("data5/" + name + "/render");
mkdir("data5/" + name + "/visual");
for i = 400 : 6 : 694
    coef = readmatrix("coef1Lime/coef" + num2str(i) + ".txt");
    id = fopen("data5/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = mallard(40.0, 4.0, 0.03);
    writematrix(xyvals, "data5/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data5/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data5/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(1, xyvals, info);

% Color 2
name = "mallard2Green";
mkdir("data5/" + name);
mkdir("data5/" + name + "/coefs");
mkdir("data5/" + name + "/geometry");
mkdir("data5/" + name + "/render");
mkdir("data5/" + name + "/visual");
for i = 400 : 6 : 694
    coef = readmatrix("coef2Green/coef" + num2str(i) + ".txt");
    id = fopen("data5/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = mallard(40.0, 4.0, 0.03);
    writematrix(xyvals, "data5/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data5/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data5/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(2, xyvals, info);

% Color 3
name = "mallard3Cyan";
mkdir("data5/" + name);
mkdir("data5/" + name + "/coefs");
mkdir("data5/" + name + "/geometry");
mkdir("data5/" + name + "/render");
mkdir("data5/" + name + "/visual");
for i = 400 : 6 : 694
    coef = readmatrix("coef3Cyan/coef" + num2str(i) + ".txt");
    id = fopen("data5/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = mallard(40.0, 4.0, 0.03);
    writematrix(xyvals, "data5/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data5/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data5/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(3, xyvals, info);

% Color 4
name = "mallard4Blue";
mkdir("data5/" + name);
mkdir("data5/" + name + "/coefs");
mkdir("data5/" + name + "/geometry");
mkdir("data5/" + name + "/render");
mkdir("data5/" + name + "/visual");
for i = 400 : 6 : 694
    coef = readmatrix("coef4Blue/coef" + num2str(i) + ".txt");
    id = fopen("data5/" + name + "/coefs/wvl" + num2str(i) + ".binary", 'w');
    fwrite(id, [size(coef, 1) size(coef, 2); 0 0], 'int');
    fwrite(id, coef, 'double');
    fclose(id);
end
for n = 1 : 50
    [xyvals, info, regularity] = mallard(40.0, 4.0, 0.03);
    writematrix(xyvals, "data5/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data5/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(regularity, "data5/" + name + "/geometry/regularity" + num2str(n) + ".txt");
end
visualize(4, xyvals, info);