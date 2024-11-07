mkdir("data4");
name = "hummingPink";
mkdir("data4/" + name);
mkdir("data4/" + name + "/geometry");
mkdir("data4/" + name + "/render");
for n = 1 : 50
    [xyvals, info] = hummingbird(40.0, 0.010, 0.090, 0.160, 0.100, 0.060, 0.020, 2.0, 12, 0.03);
    writematrix(xyvals, "data4/" + name + "/geometry/xyvals" + num2str(n) + ".txt");
    writematrix(info, "data4/" + name + "/geometry/info" + num2str(n) + ".txt");
end
visualize(1, xyvals, info);