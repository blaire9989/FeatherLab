mkdir("data2");

% Color 1
name = "starling1Aqua";
mkdir("data2/" + name);
mkdir("data2/" + name + "/geometry");
mkdir("data2/" + name + "/render");
for n = 1 : 50
    [xyvals, info, mRatio] = starling(25.0, 4.0, 0.300, 0.250, 0.03);
    writematrix(xyvals, "data2/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data2/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data2/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(1, xyvals, info);

% Color 2
name = "starling2Green";
mkdir("data2/" + name);
mkdir("data2/" + name + "/geometry");
mkdir("data2/" + name + "/render");
for n = 1 : 50
    [xyvals, info, mRatio] = starling(25.0, 4.0, 0.320, 0.250, 0.03);
    writematrix(xyvals, "data2/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data2/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data2/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(2, xyvals, info);

% Color 3
name = "starling3Yellow";
mkdir("data2/" + name);
mkdir("data2/" + name + "/geometry");
mkdir("data2/" + name + "/render");
for n = 1 : 50
    [xyvals, info, mRatio] = starling(25.0, 4.0, 0.340, 0.250, 0.03);
    writematrix(xyvals, "data2/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data2/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data2/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(3, xyvals, info);

% Color 4
name = "starling4Peach";
mkdir("data2/" + name);
mkdir("data2/" + name + "/geometry");
mkdir("data2/" + name + "/render");
for n = 1 : 50
    [xyvals, info, mRatio] = starling(25.0, 4.0, 0.360, 0.250, 0.03);
    writematrix(xyvals, "data2/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data2/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data2/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(4, xyvals, info);

% Color 5
name = "starling5Pink";
mkdir("data2/" + name);
mkdir("data2/" + name + "/geometry");
mkdir("data2/" + name + "/render");
for n = 1 : 50
    [xyvals, info, mRatio] = starling(25.0, 4.0, 0.380, 0.250, 0.03);
    writematrix(xyvals, "data2/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data2/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data2/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(5, xyvals, info);

% Color 6
name = "starling6Purple";
mkdir("data2/" + name);
mkdir("data2/" + name + "/geometry");
mkdir("data2/" + name + "/render");
for n = 1 : 50
    [xyvals, info, mRatio] = starling(25.0, 4.0, 0.400, 0.250, 0.03);
    writematrix(xyvals, "data2/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data2/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data2/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(6, xyvals, info);

% Color 7
name = "starling7Lavender";
mkdir("data2/" + name);
mkdir("data2/" + name + "/geometry");
mkdir("data2/" + name + "/render");
for n = 1 : 50
    [xyvals, info, mRatio] = starling(25.0, 4.0, 0.420, 0.250, 0.03);
    writematrix(xyvals, "data2/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data2/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data2/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(7, xyvals, info);

% Color 8
name = "starling8Teal";
mkdir("data2/" + name);
mkdir("data2/" + name + "/geometry");
mkdir("data2/" + name + "/render");
for n = 1 : 50
    [xyvals, info, mRatio] = starling(25.0, 4.0, 0.440, 0.250, 0.03);
    writematrix(xyvals, "data2/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data2/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data2/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(8, xyvals, info);