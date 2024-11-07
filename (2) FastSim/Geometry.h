#pragma once
#include "constants.h"
#include "parallel.h"

class Geometry {
    public:
        int type, nTheta, nPhi, nComponents, nPoints;
        double lambda, air, k, period;
        dcomp keratin, melanin;
        VectorXd xtop, ytop;
        MatrixXd xyvals, layer_height;
        MatrixXi info;
        MatrixXcd layer_ior;
        vector<vector<MatrixXcd>> photonic_coef;
        
        Geometry(int type, double lambda, MatrixXd photonics);
        void simulateBarbule(MatrixXd xyvals, MatrixXi info, MatrixXd melaninRatio);
        void simpleLayerStructures(MatrixXd melaninRatio);
        void hummingbirdStructures();
        Vector2d intersectLayer(double a, double b, double c, VectorXd xvec, VectorXd yvec);
        void photonicReference(MatrixXd coef);
        bool isEmptySim(MatrixXcd row_data);
};