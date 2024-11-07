#pragma once
#include "constants.h"
#include "parallel.h"
#include "Geometry.h"

class WaveSim {
    public:
        int type, nTheta, order;
        double lambda, air, k, omega, eps;
        VectorXd pvals, wvals, xtop, ytop;
        float inci_p1, inci_p2;
        VectorXcf Jt1, Jz1, Mt1, Mz1, Jt2, Jz2, Mt2, Mz2;
        MatrixXf farfield1, farfield2;
        Geometry* geom;
        
        WaveSim(Geometry* geom);
        void simulate(double theta, double phi, int resolution1, int resolution2, MatrixXd regularity);
        void computeIncidentPower(double theta, double phi);
        void computeSurfaceCurrents(double theta, double phi, MatrixXd regularity);
        Vector2cd compositeR(VectorXd height, VectorXcd ior, double tilt);
        void queryPhotonics(double theta, double phi, Vector2d pt, Vector2d tvec, Vector3cd& Er1, Vector3cd& Hr1, Vector3cd& Er2, Vector3cd& Hr2);
        void computeFarField(double theta, int resolution, bool use_projected, MatrixXf& farfield);
        MatrixXcf evaluateScattering(VectorXd xvals, VectorXd yvals, double theta, double phi, VectorXcf Jt, VectorXcf Jz, VectorXcf Mt, VectorXcf Mz);
};