mkdir("data3");

% Color 1
name = "bronzewing1Pink";
mkdir("data3/" + name);
mkdir("data3/" + name + "/geometry");
mkdir("data3/" + name + "/render");
mkdir("data3/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = bronzewing(40.0, 0.012, 0.099, 0.100, 6, 0.03);
    writematrix(xyvals, "data3/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data3/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data3/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(1, xyvals, info);

% Color 2
name = "bronzewing2Red";
mkdir("data3/" + name);
mkdir("data3/" + name + "/geometry");
mkdir("data3/" + name + "/render");
mkdir("data3/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = bronzewing(40.0, 0.012, 0.091, 0.098, 6, 0.03);
    writematrix(xyvals, "data3/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data3/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data3/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(2, xyvals, info);

% Color 3
name = "bronzewing3Orange";
mkdir("data3/" + name);
mkdir("data3/" + name + "/geometry");
mkdir("data3/" + name + "/render");
mkdir("data3/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = bronzewing(40.0, 0.012, 0.089, 0.091, 6, 0.03);
    writematrix(xyvals, "data3/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data3/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data3/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(3, xyvals, info);

% Color 4
name = "bronzewing4Yellow";
mkdir("data3/" + name);
mkdir("data3/" + name + "/geometry");
mkdir("data3/" + name + "/render");
mkdir("data3/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = bronzewing(40.0, 0.012, 0.084, 0.087, 6, 0.03);
    writematrix(xyvals, "data3/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data3/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data3/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(4, xyvals, info);

% Color 5
name = "bronzewing5Green";
mkdir("data3/" + name);
mkdir("data3/" + name + "/geometry");
mkdir("data3/" + name + "/render");
mkdir("data3/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = bronzewing(40.0, 0.012, 0.082, 0.077, 6, 0.03);
    writematrix(xyvals, "data3/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data3/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data3/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(5, xyvals, info);

% Color 6
name = "bronzewing6Cyan";
mkdir("data3/" + name);
mkdir("data3/" + name + "/geometry");
mkdir("data3/" + name + "/render");
mkdir("data3/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = bronzewing(40.0, 0.012, 0.077, 0.072, 6, 0.03);
    writematrix(xyvals, "data3/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data3/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data3/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(6, xyvals, info);

% Color 7
name = "bronzewing7Blue";
mkdir("data3/" + name);
mkdir("data3/" + name + "/geometry");
mkdir("data3/" + name + "/render");
mkdir("data3/" + name + "/visual");
for n = 1 : 50
    [xyvals, info, mRatio] = bronzewing(40.0, 0.012, 0.075, 0.063, 6, 0.03);
    writematrix(xyvals, "data3/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data3/" + name + "/geometry/info" + num2str(n) + ".txt");
    writematrix(mRatio, "data3/" + name + "/geometry/melaninRatio" + num2str(n) + ".txt");
end
visualize(7, xyvals, info);