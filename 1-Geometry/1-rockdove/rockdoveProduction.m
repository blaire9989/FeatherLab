mkdir("data1");

% Color 1
name = "rockdove1Forest";
mkdir("data1/" + name);
mkdir("data1/" + name + "/geometry");
mkdir("data1/" + name + "/render");
mkdir("data1/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = rockdove(20.0, 4.0, 0.595, 0.03);
    writematrix(xyvals, "data1/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data1/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data1/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(1, xyvals, info);

% Color 2
name = "rockdove2Jade";
mkdir("data1/" + name);
mkdir("data1/" + name + "/geometry");
mkdir("data1/" + name + "/render");
mkdir("data1/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = rockdove(20.0, 4.0, 0.582, 0.03);
    writematrix(xyvals, "data1/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data1/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data1/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(2, xyvals, info);

% Color 3
name = "rockdove3Sea";
mkdir("data1/" + name);
mkdir("data1/" + name + "/geometry");
mkdir("data1/" + name + "/render");
mkdir("data1/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = rockdove(20.0, 4.0, 0.569, 0.03);
    writematrix(xyvals, "data1/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data1/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data1/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(3, xyvals, info);

% Color 4
name = "rockdove4Teal";
mkdir("data1/" + name);
mkdir("data1/" + name + "/geometry");
mkdir("data1/" + name + "/render");
mkdir("data1/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = rockdove(20.0, 4.0, 0.556, 0.03);
    writematrix(xyvals, "data1/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data1/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data1/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(4, xyvals, info);

% Color 5
name = "rockdove5Violet";
mkdir("data1/" + name);
mkdir("data1/" + name + "/geometry");
mkdir("data1/" + name + "/render");
mkdir("data1/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = rockdove(20.0, 4.0, 0.543, 0.03);
    writematrix(xyvals, "data1/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data1/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data1/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(5, xyvals, info);

% Color 6
name = "rockdove6Purple";
mkdir("data1/" + name);
mkdir("data1/" + name + "/geometry");
mkdir("data1/" + name + "/render");
mkdir("data1/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = rockdove(20.0, 4.0, 0.530, 0.03);
    writematrix(xyvals, "data1/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data1/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data1/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(6, xyvals, info);