#include <unistd.h>
#include "BRDF.h"
#include "WaveSim.h"

int main(int argc, char **argv) {
    int opt, iTheta, iPhi, type, nInstance, visualize;
    double lambda;
    string xyname;
    while ((opt = getopt(argc, argv, "a:b:i:l:n:v:z:")) != -1) {
        switch (opt) {
            case 'a': iTheta = atoi(optarg); break;
            case 'b': iPhi = atoi(optarg); break;
            case 'i': type = atoi(optarg); break;
            case 'l': lambda = atof(optarg) * 0.001; break;
            case 'n': nInstance = atoi(optarg); break;
            case 'v': visualize = atoi(optarg); break;
            case 'z': xyname = optarg; break;
        }
    }
    int dim = 20, resolution1 = 400, resolution2 = 1024;
    double theta = acos(1.0 / dim * iTheta);
    double phi = acos(1.0 / dim * iPhi);
    MatrixXd photonics, melaninRatio, regularity;
    if (type >= 5)
        readBinary("data" + to_string(type) + "/" + xyname + "/coefs/wvl" + to_string((int)round(1000 * lambda)) + ".binary", photonics);
    Geometry* geom = new Geometry(type, lambda, photonics);
    WaveSim* sim = new WaveSim(geom);
    BRDF* brdf = new BRDF(resolution1, resolution2, nInstance);
    
    // Loop through all barbule instances
    VectorXd offset = VectorXd::LinSpaced(nInstance, -2.0 * M_PI / 180.0, 2.0 * M_PI / 180.0);
    for (int i = 1; i <= nInstance; i++) {
        MatrixXd xyvals = readData("data" + to_string(type) + "/" + xyname + "/geometry/xyvals" + to_string(i) + ".txt");
        MatrixXi info = readData("data" + to_string(type) + "/" + xyname + "/geometry/info" + to_string(i) + ".txt").cast<int>();
        if (type <= 3)
            melaninRatio = readData("data" + to_string(type) + "/" + xyname + "/geometry/melaninRatio" + to_string(i) + ".txt");
        if (type >= 5)
            regularity = readData("data" + to_string(type) + "/" + xyname + "/geometry/regularity" + to_string(i) + ".txt");
        geom->simulateBarbule(xyvals, info, melaninRatio);
        sim->simulate(theta, phi + offset(i - 1), resolution1, resolution2, regularity);
        brdf->data1.block(0, 14 * (i - 1), resolution1, 14) = sim->farfield1;
        brdf->data2.block(0, 14 * (i - 1), resolution2, 14) = sim->farfield2;
    }

    // Compute the mean BRDF parameters and field value statistics
    int index = dim * (iTheta - 1) + iPhi;
    brdf->computeParameters(theta, phi);
    if (visualize == 0) {
        VectorXd parameters = VectorXd::Zero(12);
        parameters(0) = brdf->m1;
        parameters(1) = brdf->phi1;
        parameters.block(2, 0, 3, 1) = brdf->f1;
        parameters(5) = brdf->m2;
        parameters(6) = brdf->phi2;
        parameters.block(7, 0, 3, 1) = brdf->f2;
        parameters(10) = brdf->rate;
        parameters(11) = brdf->l;
        writeData("data" + to_string(type) + "/" + xyname + "/render/param_" + to_string((int)round(1000 * lambda)) + "_" + to_string(index) + ".binary", parameters);
        if (index == dim * dim && lambda < 0.405)
            writeData("data" + to_string(type) + "/" + xyname + "/render/noise.binary", brdf->noise);
    } else {
        MatrixXf visual = MatrixXf::Zero(resolution1, 5);
        visual.block(0, 0, resolution1, 1) = brdf->brdf_sim;
        visual.block(0, 1, resolution1, 1) = brdf->brdf_fit;
        visual.block(0, 2, resolution1, 1) = 0.5f * brdf->data1.block(0, 0, resolution1, 1) + 0.5f * brdf->data1.block(0, 7, resolution1, 1);
        visual.block(0, 3, resolution1, 1) = 0.5f * brdf->data1.block(0, 14, resolution1, 1) + 0.5f * brdf->data1.block(0, 21, resolution1, 1);
        visual.block(0, 4, resolution1, 1) = 0.5f * brdf->data1.block(0, 28, resolution1, 1) + 0.5f * brdf->data1.block(0, 35, resolution1, 1);
        writeData("data" + to_string(type) + "/" + xyname + "/visual/pattern_" + to_string((int)round(1000 * lambda)) + "_" + to_string(index) + ".binary", visual);
    }
    delete geom;
    delete sim;
    delete brdf;
}