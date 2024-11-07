/* This module computes all the BRDF parameters describing our reflectance models for iridescent feathers.
Users are welcome to read but do not need to fully understand methods in this file.
PLEASE DO NOT MODIFY THE FOLLOWING CODE. */

#include "BRDF.h"

/// @brief Allocating data structures for storing scattering quantities computed from each barbule instance
/// @param resolution1, resolution2: specified resolutions used for computing the scattering quantities
/// @param nInstance: the number of barbule instances we simulate
BRDF::BRDF(int resolution1, int resolution2, int nInstance) {
    this->resolution1 = resolution1;
    this->resolution2 = resolution2;
    this->nInstance = nInstance;
    data1 = MatrixXf::Zero(resolution1, 14 * nInstance);
    data2 = MatrixXf::Zero(resolution2, 14 * nInstance);
}

/// @brief The umbrella function for computing all the required BRDF parameters
/// @param theta_i, phi_i: specifying the direction of the incident light
void BRDF::computeParameters(double theta_i, double phi_i) {
    brdf_sim = VectorXf::Zero(resolution1);
    for (int i = 0; i < nInstance; i++) {
        brdf_sim = brdf_sim + 0.5 * data1.block(0, 14 * i, resolution1, 1) / nInstance;
        brdf_sim = brdf_sim + 0.5 * data1.block(0, 14 * i + 7, resolution1, 1) / nInstance;
    }
    fitGGXParameters(theta_i, phi_i, brdf_sim, brdf_fit, m1, phi1);
    signal_sim = 0.5 * computeSignalLevel(1) + 0.5 * computeSignalLevel(2);
    fitGGXParameters(theta_i, phi_i, signal_sim, signal_fit, m2, phi2);
    f1 = VectorXd::Zero(3);
    for (int i = 0; i < 3; i++)
        f1(i) = computeBRDFFresnel(theta_i, m1, phi1, 0.02 * (i + 2)) * cos(theta_i) * cos(phi_i);
    double s0 = signal_sim.sum() / brdf_sim.sum();
    f2 = VectorXd::Zero(3);
    for (int i = 0; i < 3; i++)
        f2(i) = computeSignalFresnel(theta_i, 0.02 * (i + 2), f1(i), s0);
    computeEnergyRates();
    noise = MatrixXf::Zero(resolution2, 4 * nInstance);
    generateNoiseData(1);
    generateNoiseData(2);
    l = correlationLength(theta_i, m1, phi1);
}

/// @brief Fit the average reflectance function to a GGX model
/// @param theta_i, phi_i: specifying the direction of the incident light
/// @param original, fitted: the original and fitted average reflectance function
/// @param m, phi_c: the GGX model parameters, see Section 5.2 of our paper
void BRDF::fitGGXParameters(double theta_i, double phi_i, VectorXf original, VectorXf& fitted, double& m, double& phi_c) {
    double value = cos(theta_i) * cos(phi_i);
    if (value <= 0.05)
        phi_c = phi_i;
    else if (value >= 0.10)
        phi_c = findPhi(theta_i, phi_i, original);
    else {
        double s = (value - 0.05) / 0.05;
        phi_c = findPhi(theta_i, phi_i, original);
        phi_c = (1 - s) * phi_i + s * phi_c;
    }
    m = findM(theta_i, phi_c, original, fitted);
}

/// @brief Optimize the phi_eff parameter in the GGX model fitting
/// @param theta_i, phi_i: specifying the direction of the incident light
/// @param original: the original average reflectance function that we hope to fit
double BRDF::findPhi(double theta_i, double phi_i, VectorXf original) {
    VectorXd mvec1 = VectorXd::LinSpaced(20, 0.025, 0.975);
    VectorXd mvec2 = VectorXd::LinSpaced(200, 0.005, 1.000);
    VectorXd mvec3 = VectorXd::LinSpaced(1000, 0.001, 1.000);
    int nPhi = 180;
    VectorXd phivec = VectorXd::LinSpaced(nPhi, -0.5 * M_PI + 0.5 * M_PI / nPhi, 0.5 * M_PI - 0.5 * M_PI / nPhi);
    int ind1, ind2;
    double m = -1.0;
    float energy = original.sum();
    double phi_c;

    // Round 1 fitting
    float minimum = 1e6;
    for (int i = 0; i < mvec1.size(); i++) {
        VectorXf residual = VectorXf::Zero(nPhi);
        parallel_for(nPhi, [&](int start, int end) {
        for (int j = start; j < end; j++) {
            VectorXf G = GGX(mvec1(i), theta_i, phivec(j), resolution1);
            G = G / G.sum() * energy;
            residual(j) = (G - original).cwiseAbs2().sum();
        }
        } );
        if (residual.minCoeff() < minimum) {
            minimum = residual.minCoeff();
            m = mvec1(i);
        }
    }

    // Round 2 fitting
    refineParameterRange(ind1, ind2, mvec2, m, 0.05);
    minimum = 1e6;
    for (int i = ind1; i <= ind2; i++) {
        VectorXf residual = VectorXf::Zero(nPhi);
        parallel_for(nPhi, [&](int start, int end) {
        for (int j = start; j < end; j++) {
            VectorXf G = GGX(mvec2(i), theta_i, phivec(j), resolution1);
            G = G / G.sum() * energy;
            residual(j) = (G - original).cwiseAbs2().sum();
        }
        } );
        if (residual.minCoeff() < minimum) {
            minimum = residual.minCoeff();
            m = mvec2(i);
        }
    }

    // Round 3 fitting
    refineParameterRange(ind1, ind2, mvec3, m, 0.005);
    minimum = 1e6;
    for (int i = ind1; i <= ind2; i++) {
        VectorXf residual = VectorXf::Zero(nPhi);
        parallel_for(nPhi, [&](int start, int end) {
        for (int j = start; j < end; j++) {
            VectorXf G = GGX(mvec3(i), theta_i, phivec(j), resolution1);
            G = G / G.sum() * energy;
            residual(j) = (G - original).cwiseAbs2().sum();
        }
        } );
        for (int j = 0; j < nPhi; j++) {
            if (residual(j) < minimum) {
                minimum = residual(j);
                phi_c = phivec(j);
            }
        }
    }
    return phi_c;
}

/// @brief Optimize the m_eff parameter in the GGX model fitting
/// @param theta_i: specifying the longitudinal angle of the incident light
/// @param phi_c: the chosen phi_eff parameter in the GGX model
/// @param original, fitted: the original and fitted average reflectance function
double BRDF::findM(double theta_i, double phi_c, VectorXf original, VectorXf& fitted) {
    VectorXd mvec1 = VectorXd::LinSpaced(20, 0.025, 0.975);
    VectorXd mvec2 = VectorXd::LinSpaced(200, 0.005, 1.000);
    VectorXd mvec3 = VectorXd::LinSpaced(1000, 0.001, 1.000);
    int ind1, ind2;
    double m = -1.0;
    float energy = original.sum();

    // Round 1 fitting
    float minimum = 1e6;
    for (int i = 0; i < mvec1.size(); i++) {
        VectorXf G = GGX(mvec1(i), theta_i, phi_c, resolution1);
        G = G / G.sum() * energy;
        float residual = (G - original).cwiseAbs2().sum();
        if (residual < minimum) {
            minimum = residual;
            m = mvec1(i);
        }
    }

    // Round 2 fitting
    refineParameterRange(ind1, ind2, mvec2, m, 0.05);
    minimum = 1e6;
    for (int i = ind1; i <= ind2; i++) {
        VectorXf G = GGX(mvec2(i), theta_i, phi_c, resolution1);
        G = G / G.sum() * energy;
        float residual = (G - original).cwiseAbs2().sum();
        if (residual < minimum) {
            minimum = residual;
            m = mvec2(i);
        }
    }

    // Round 3 fitting
    refineParameterRange(ind1, ind2, mvec3, m, 0.005);
    minimum = 1e6;
    for (int i = ind1; i <= ind2; i++) {
        VectorXf G = GGX(mvec3(i), theta_i, phi_c, resolution1);
        G = G / G.sum() * energy;
        float residual = (G - original).cwiseAbs2().sum();
        if (residual < minimum) {
            minimum = residual;
            m = mvec3(i);
        }
    }
    VectorXf G = GGX(m, theta_i, phi_c, resolution1);
    fitted = G / G.sum() * energy;
    return m;
}

/// @brief Evaluate the GGX model given the model parameters
/// @param m, theta_i, phi_i: the roughness parameter and the incident direction
/// @param resolution: the resolution at which we evaluate the GGX model
VectorXf BRDF::GGX(double m, double theta_i, double phi_i, int resolution) {
    Vector3d wi(cos(theta_i) * sin(phi_i), sin(theta_i), cos(theta_i) * cos(phi_i));
    VectorXd phi_o = VectorXd::LinSpaced(resolution, -0.5 * M_PI + 0.5 * M_PI / resolution, 0.5 * M_PI - 0.5 * M_PI / resolution);
    VectorXf G = VectorXf::Zero(resolution);
    for (int i = 0; i < resolution; i++) {
        Vector3d wo(cos(theta_i) * sin(phi_o(i)), -sin(theta_i), cos(theta_i) * cos(phi_o(i)));
        Vector3d h = (wi + wo) / (wi + wo).norm();
        double D = pow(m / (m * m + (1 - m * m) * h(0) * h(0)), 2) / M_PI; 
        double Gi = 2.0 / (1.0 + sqrt(1 + m * m * (1 - wi(2) * wi(2)) / (wi(2) * wi(2))));
        double Go = 2.0 / (1.0 + sqrt(1 + m * m * (1 - wo(2) * wo(2)) / (wo(2) * wo(2))));
        G(i) = (float)(D * Gi * Go);
    }
    return G;
}

/// @brief Refine the range of m_eff parameters for more accurate fitting
/// @param ind_min, ind_max, range: specifying the narrowed range of m_eff
/// @param val, delta: specifying the current, rough estimate of m_eff
void BRDF::refineParameterRange(int& ind_min, int& ind_max, VectorXd range, double val, double delta) {
    ind_min = -1;
    ind_max = -1;
    double minimum1 = 1e6, minimum2 = 1e6;
    for (int i = 0; i < range.size(); i++) {
        double diff1 = abs(range(i) - (val - delta));
        if (diff1 < minimum1) {
            minimum1 = diff1;
            ind_min = i;
        }
        double diff2 = abs(range(i) - (val + delta));
        if (diff2 <= minimum2) {
            minimum2 = diff2;
            ind_max = i;
        }
    }
}

/// @brief Compute the statistics function S described in Section 5.3 of our paper
/// @param polarization: indicates which polarization we are computing the function S for (1 or 2). S is computed for both polarizations and averaged
VectorXf BRDF::computeSignalLevel(int polarization) {
    int offset = 7 * (polarization - 1);
    VectorXf mu_xR = VectorXf::Zero(resolution1);
    VectorXf mu_xI = VectorXf::Zero(resolution1);
    VectorXf mu_yR = VectorXf::Zero(resolution1);
    VectorXf mu_yI = VectorXf::Zero(resolution1);
    VectorXf mu_zR = VectorXf::Zero(resolution1);
    VectorXf mu_zI = VectorXf::Zero(resolution1);

    // Compute mean field values
    for (int i = 0; i < nInstance; i++) {
        mu_xR = mu_xR + data1.block(0, 14 * i + 1 + offset, resolution1, 1) / nInstance;
        mu_xI = mu_xI + data1.block(0, 14 * i + 2 + offset, resolution1, 1) / nInstance;
        mu_yR = mu_yR + data1.block(0, 14 * i + 3 + offset, resolution1, 1) / nInstance;
        mu_yI = mu_yI + data1.block(0, 14 * i + 4 + offset, resolution1, 1) / nInstance;
        mu_zR = mu_zR + data1.block(0, 14 * i + 5 + offset, resolution1, 1) / nInstance;
        mu_zI = mu_zI + data1.block(0, 14 * i + 6 + offset, resolution1, 1) / nInstance;
    }

    // Compute the sum of squares of mean field values
    double eps = 1 / (mu * c0 * c0);
    float scale = 0.5 * sqrt(eps / mu);
    VectorXf signal = mu_xR.cwiseAbs2() + mu_xI.cwiseAbs2() + mu_yR.cwiseAbs2() + mu_yI.cwiseAbs2() + mu_zR.cwiseAbs2() + mu_zI.cwiseAbs2();
    return scale * signal;
}

/// @brief Compute the effective Fresnel parameter f_eff in Section 5.2 of our paper
/// @param theta_i: specifying the longitudinal angle of the incident light
/// @param m, phi_c: the chosen m_eff and phi_eff parameters
/// @param sigma: the chosen sigma parameter for the y-axis spread in our full hemisphere BRDF model
double BRDF::computeBRDFFresnel(double theta_i, double m, double phi_c, double sigma) {
    int resAccurate = 2000;
    double energy = cos(theta_i) * brdf_sim.sum() * M_PI / resolution1;
    Vector3d wi(cos(theta_i) * sin(phi_c), sin(theta_i), cos(theta_i) * cos(phi_c));
    double Gi = 2.0 / (1.0 + sqrt(1 + m * m * (1 - wi(2) * wi(2)) / (wi(2) * wi(2))));

    // Compute the effective Fresnel term: a scale factor
    double delta = M_PI / resAccurate;
    VectorXd theta_o = VectorXd::LinSpaced(resAccurate, -0.5 * M_PI + 0.5 * delta, 0.5 * M_PI - 0.5 * delta);
    VectorXd phi_o = VectorXd::LinSpaced(resAccurate, -0.5 * M_PI + 0.5 * delta, 0.5 * M_PI - 0.5 * delta);
    double total = 0;
    for (int i = 0; i < resAccurate; i++) {
        double scale = delta * delta * cos(theta_o(i));
        for (int j = 0; j < resAccurate; j++) {
            Vector3d wo(cos(theta_o(i)) * sin(phi_o(j)), sin(theta_o(i)), cos(theta_o(i)) * cos(phi_o(j)));
            double Go = 2.0 / (1.0 + sqrt(1 + m * m * (1 - wo(2) * wo(2)) / (wo(2) * wo(2))));
            Vector3d h = (wi + wo) / (wi + wo).norm();
            double D1 = pow(m / (m * m + (1 - m * m) * h(0) * h(0)), 2) / M_PI;
            double D2 = exp(-h(1) * h(1) / (sigma * sigma));
            total = total + scale * D1 * D2 * Gi * Go / (4.0 * wi(2));
        }
    }
    return energy / total;
}

/// @brief Compute the effective Fresnel parameter f_sta in Section 5.3 of our paper
/// @param theta_i: specifying the longitudinal angle of the incident light
/// @param sigma: the chosen sigma parameter for the y-axis spread in our full hemisphere BRDF model
/// @param f: the computed f_eff parameter
/// @param s0: a signal-to-noise estimate that characterizes the variances of reflectance distributions among barbule instances
double BRDF::computeSignalFresnel(double theta_i, double sigma, double f, double s0) {
    int resAccurate = 2000;
    Vector3d wi1(cos(theta_i) * sin(phi1), sin(theta_i), cos(theta_i) * cos(phi1));
    Vector3d wi2(cos(theta_i) * sin(phi2), sin(theta_i), cos(theta_i) * cos(phi2));
    double Gi1 = 2.0 / (1.0 + sqrt(1 + m1 * m1 * (1 - wi1(2) * wi1(2)) / (wi1(2) * wi1(2))));
    double Gi2 = 2.0 / (1.0 + sqrt(1 + m2 * m2 * (1 - wi2(2) * wi2(2)) / (wi2(2) * wi2(2))));

    // Compute the scale factor on the signal level
    VectorXd phi_o = VectorXd::LinSpaced(resAccurate, -0.5 * M_PI + 0.5 * M_PI / resAccurate, 0.5 * M_PI - 0.5 * M_PI / resAccurate);
    double total1 = 0, total2 = 0;
    for (int i = 0; i < resAccurate; i++) {
        Vector3d wo(cos(theta_i) * sin(phi_o(i)), -sin(theta_i), cos(theta_i) * cos(phi_o(i)));
        double Go1 = 2.0 / (1.0 + sqrt(1 + m1 * m1 * (1 - wo(2) * wo(2)) / (wo(2) * wo(2))));
        double Go2 = 2.0 / (1.0 + sqrt(1 + m2 * m2 * (1 - wo(2) * wo(2)) / (wo(2) * wo(2))));
        Vector3d h1 = (wi1 + wo) / (wi1 + wo).norm();
        Vector3d h2 = (wi2 + wo) / (wi2 + wo).norm();
        double D11 = pow(m1 / (m1 * m1 + (1 - m1 * m1) * h1(0) * h1(0)), 2) / M_PI;
        double D21 = exp(-h1(1) * h1(1) / (sigma * sigma));
        total1 = total1 + f * D11 * D21 * Gi1 * Go1 / (4.0 * wi1(2));
        double D12 = pow(m2 / (m2 * m2 + (1 - m2 * m2) * h2(0) * h2(0)), 2) / M_PI;
        double D22 = exp(-h2(1) * h2(1) / (sigma * sigma));
        total2 = total2 + D12 * D22 * Gi2 * Go2 / (4.0 * wi2(2));
    }
    return s0 * total1 / total2;
}

/// @brief Compute a parameter r not discussed in our paper, used for distinguishing between the two polarizations when synthesizing new BRDF instances
void BRDF::computeEnergyRates() {
    VectorXf intensity1 = VectorXf::Zero(resolution1);
    VectorXf intensity2 = VectorXf::Zero(resolution1);
    for (int i = 0; i < nInstance; i++) {
        intensity1 = intensity1 + data1.block(0, 14 * i, resolution1, 1);
        intensity2 = intensity2 + data1.block(0, 14 * i + 7, resolution1, 1);
    }
    rate = intensity1.sum() / (intensity1.sum() + intensity2.sum());
}

/// @brief Generate some noise data for computing more noise functions used in synthesized BRDF instances
/// @param polarization: indicates which polarization we are computing the noise data from. Both polarizations are used
void BRDF::generateNoiseData(int polarization) {
    int offset = 7 * (polarization - 1);
    VectorXf mu_xR = VectorXf::Zero(resolution2);
    VectorXf mu_xI = VectorXf::Zero(resolution2);
    VectorXf mu_yR = VectorXf::Zero(resolution2);
    VectorXf mu_yI = VectorXf::Zero(resolution2);
    VectorXf mu_zR = VectorXf::Zero(resolution2);
    VectorXf mu_zI = VectorXf::Zero(resolution2);
    VectorXf var_xR = VectorXf::Zero(resolution2);
    VectorXf var_xI = VectorXf::Zero(resolution2);
    VectorXf var_yR = VectorXf::Zero(resolution2);
    VectorXf var_yI = VectorXf::Zero(resolution2);
    VectorXf var_zR = VectorXf::Zero(resolution2);
    VectorXf var_zI = VectorXf::Zero(resolution2);

    // Compute the mean and variance of field values
    for (int i = 0; i < nInstance; i++) {
        mu_xR = mu_xR + data2.block(0, 14 * i + 1 + offset, resolution2, 1) / nInstance;
        mu_xI = mu_xI + data2.block(0, 14 * i + 2 + offset, resolution2, 1) / nInstance;
        mu_yR = mu_yR + data2.block(0, 14 * i + 3 + offset, resolution2, 1) / nInstance;
        mu_yI = mu_yI + data2.block(0, 14 * i + 4 + offset, resolution2, 1) / nInstance;
        mu_zR = mu_zR + data2.block(0, 14 * i + 5 + offset, resolution2, 1) / nInstance;
        mu_zI = mu_zI + data2.block(0, 14 * i + 6 + offset, resolution2, 1) / nInstance;
    }
    for (int i = 0; i < nInstance; i++) {
        var_xR = var_xR + (data2.block(0, 14 * i + 1 + offset, resolution2, 1) - mu_xR).cwiseAbs2() / (nInstance - 1);
        var_xI = var_xI + (data2.block(0, 14 * i + 2 + offset, resolution2, 1) - mu_xI).cwiseAbs2() / (nInstance - 1);
        var_yR = var_yR + (data2.block(0, 14 * i + 3 + offset, resolution2, 1) - mu_yR).cwiseAbs2() / (nInstance - 1);
        var_yI = var_yI + (data2.block(0, 14 * i + 4 + offset, resolution2, 1) - mu_yI).cwiseAbs2() / (nInstance - 1);
        var_zR = var_zR + (data2.block(0, 14 * i + 5 + offset, resolution2, 1) - mu_zR).cwiseAbs2() / (nInstance - 1);
        var_zI = var_zI + (data2.block(0, 14 * i + 6 + offset, resolution2, 1) - mu_zI).cwiseAbs2() / (nInstance - 1);
    }

    // Select the largest field value component among x, y, z and compute the normalized noise data
    float energyX = (mu_xR.cwiseAbs2() + mu_xI.cwiseAbs2() + var_xR + var_xI).sum();
    float energyY = (mu_yR.cwiseAbs2() + mu_yI.cwiseAbs2() + var_yR + var_yI).sum();
    float energyZ = (mu_zR.cwiseAbs2() + mu_zI.cwiseAbs2() + var_zR + var_zI).sum();
    if (energyX > energyY && energyX > energyZ) {
        for (int i = 0; i < resolution2; i++) {
            for (int j = 0; j < nInstance; j++) {
                noise(i, (2 * polarization - 2) * nInstance + j) = (data2(i, 14 * j + 1 + offset) - mu_xR(i)) / sqrt(var_xR(i));
                noise(i, (2 * polarization - 1) * nInstance + j) = (data2(i, 14 * j + 2 + offset) - mu_xI(i)) / sqrt(var_xI(i));
            }
        }
    } else if (energyY > energyX && energyY > energyZ) {
        for (int i = 0; i < resolution2; i++) {
            for (int j = 0; j < nInstance; j++) {
                noise(i, (2 * polarization - 2) * nInstance + j) = (data2(i, 14 * j + 3 + offset) - mu_yR(i)) / sqrt(var_yR(i));
                noise(i, (2 * polarization - 1) * nInstance + j) = (data2(i, 14 * j + 4 + offset) - mu_yI(i)) / sqrt(var_yI(i));
            }
        }
    } else {
        for (int i = 0; i < resolution2; i++) {
            for (int j = 0; j < nInstance; j++) {
                noise(i, (2 * polarization - 2) * nInstance + j) = (data2(i, 14 * j + 5 + offset) - mu_zR(i)) / sqrt(var_zR(i));
                noise(i, (2 * polarization - 1) * nInstance + j) = (data2(i, 14 * j + 6 + offset) - mu_zI(i)) / sqrt(var_zI(i));
            }
        }
    }
}

/// @brief Compute the length parameter l in Section 5.3 of our paper
/// @param theta_i: specifying the longitudinal angle of the incident light
/// @param m, phi_c: the chosen m_eff and phi_eff parameter in our BRDF model
double BRDF::correlationLength(double theta_i, double m, double phi_c) {
    Vector3d wi(cos(theta_i) * sin(phi_c), sin(theta_i), cos(theta_i) * cos(phi_c));
    double Gi = 2.0 / (1.0 + sqrt(1 + m * m * (1 - wi(2) * wi(2)) / (wi(2) * wi(2))));

    // Cut off the zero value and weight by the BRDF value
    double y0 = -sin(theta_i);
    int minInd = resolution2, maxInd = -1;
    for (int i = 0; i < resolution2; i++) {
        double x0 = -1.0 + 1.0 / resolution2 + i * 2.0 / resolution2;
        if (x0 * x0 + y0 * y0 <= 1) {
            if (i < minInd)
                minInd = i;
            if (i > maxInd)
                maxInd = i;
            Vector3d wo(x0, y0, sqrt(1 - x0 * x0 - y0 * y0));
            Vector3d h = (wi + wo) / (wi + wo).norm();
            double D = pow(m / (m * m + (1 - m * m) * h(0) * h(0)), 2) / M_PI;
            double Go = 2.0 / (1.0 + sqrt(1 + m * m * (1 - wo(2) * wo(2)) / (wo(2) * wo(2))));
            noise.block(i, 0, 1, 4 * nInstance) = D * Gi * Go * noise.block(i, 0, 1, 4 * nInstance);
        }
    }
    MatrixXf patterns = noise.block(minInd, 0, maxInd - minInd + 1, 4 * nInstance);

    // Compute the average ACF of all the noise patterns
    int N = patterns.rows();
    MatrixXf acv = MatrixXf::Zero(N / 2, 4 * nInstance);
    for (int i = 0; i < 4 * nInstance; i++) {
        float average = patterns.block(0, i, N, 1).mean();
        for (int j = 0; j < N / 2; j++) {
            float total = 0;
            for (int k = 0; k < N; k++) {
                int index = k + j;
                if (index >= N)
                    index = index - N;
                total = total + (patterns(k, i) - average) * (patterns(index, i) - average);
            }
            acv(j, i) = total / (N - 1);
        }
    }
    VectorXf mu_acv = VectorXf::Zero(N / 2);
    for (int i = 0; i < 4 * nInstance; i++)
        mu_acv = mu_acv + acv.block(0, i, N / 2, 1);
    mu_acv = mu_acv / mu_acv(0);

    // Compute autocorrelation length
    double dr = 2.0 / resolution2;
    VectorXd r = VectorXd::LinSpaced(N / 2, 0, (N / 2 - 1) * dr);
    double thres = exp(-1.0);
    int index = 0;
    while (mu_acv(index) >= thres)
        index = index + 1;
    double s = (thres - mu_acv(index - 1)) / (mu_acv(index) - mu_acv(index - 1));
    double l = (1 - s) * r(index - 1) + s * r(index);
    return l;
}