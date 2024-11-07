#pragma once
#include "constants.h"
#include "parallel.h"
#include "BRDF.h"

class BRDF {
    public:
        int resolution1, resolution2, nInstance;
        double m1, phi1, m2, phi2, rate, l;
        VectorXf brdf_sim, brdf_fit, signal_sim, signal_fit;
        VectorXd f1, f2;
        MatrixXf data1, data2, noise;
        
        BRDF(int resolution1, int resolution2, int nInstance);
        void computeParameters(double theta_i, double phi_i);
        void fitGGXParameters(double theta_i, double phi_i, VectorXf original, VectorXf& fitted, double& m, double& phi_c);
        double findPhi(double theta_i, double phi_i, VectorXf original);
        double findM(double theta_i, double phi_c, VectorXf original, VectorXf& fitted);
        VectorXf GGX(double m, double theta_i, double phi_i, int resolution);
        void refineParameterRange(int& ind_min, int& ind_max, VectorXd range, double val, double delta);
        VectorXf computeSignalLevel(int polarization);
        double computeBRDFFresnel(double theta_i, double m, double phi_c, double sigma);
        double computeSignalFresnel(double theta_i, double sigma, double f, double s0);
        void computeEnergyRates();
        void generateNoiseData(int polarization);
        double correlationLength(double theta_i, double m, double phi_c);
};