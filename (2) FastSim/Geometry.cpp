/* This module processes the specified, nano-scale geometry that describes each given barbule instance.
Layer structures in different types of barbules are analyzed, and reflection coefficients of photonic crystal structures are stored as tables.
Users are welcome to read but do not need to fully understand methods in this file.
PLEASE DO NOT MODIFY THE FOLLOWING CODE. */

#include "Geometry.h"

/// @brief Initialize the IORs of keratin and melanin. Also initialize the photonic crystal reflection coefficient table if needed
/// @param type: the type of the given barbule
/// @param type == 1: rock dove
/// @param type == 2: European starling
/// @param type == 3: bronzewing
/// @param type == 4: hummingbird
/// @param type == 5: mallard
/// @param type == 6: magpie
/// @param type == 7: peacock
/// @param lambda: the simulated wavelength
/// @param photonics: the input photonics crystal reflection coefficient table, in its raw form--it is an empty matrix if type <= 4
Geometry::Geometry(int type, double lambda, MatrixXd photonics) {
    this->type = type;
    this->lambda = lambda;
    air = 1.0;
    k = 2.0 * M_PI * air / lambda;
    keratin = 1.532 + 5890 / (lambda * lambda * 1e6);
    melanin = 1.648 + 23700 / (lambda * lambda * 1e6) + cunit * 0.56 * exp(-lambda * 1e3 / 270);
    
    // Initialize the photonic crystal reflection coefficient table if needed
    if (type >= 5) {
        photonic_coef.clear();
        this->period = photonics(0, 0);
        this->nTheta = int(photonics(0, 1));
        this->nPhi = int(photonics(0, 2)) + 1;
        MatrixXd coef = photonics.block(1, 0, 12, photonics.cols());
        photonicReference(coef);
    }
}

/// @brief Analyze the layer heights and layer IORs in the imperfect multilayer type barbules
/// @param xyvals: nanostructure of a given barbule, in the form of a 2D cross section
/// @param info: additional information, used as indices into the xy-values specifying the 2D cross section
/// @param melaninRatio: a spatially varying quantity, modeling the percentage of melanin in regions where melanosomes are dispersed in keratin or air
void Geometry::simulateBarbule(MatrixXd xyvals, MatrixXi info, MatrixXd melaninRatio) {
    this->xyvals = xyvals;
    this->info = info;
    nComponents = info.rows();
    nPoints = info(0, 1) - 2;
    xtop = xyvals.block(0, 0, info(0, 1), 1);
    ytop = xyvals.block(0, 1, info(0, 1), 1);
    layer_height = MatrixXd::Zero(nComponents - 1, nPoints);
    layer_ior = MatrixXcd::Zero(nComponents + 1, nPoints);
    
    // Slightly different simulations for different barbules
    if (type <= 3)
        simpleLayerStructures(melaninRatio);
    else if (type == 4)
        hummingbirdStructures();
}

/// @brief Analyze the spatially varying multilayer structures in rock dove, starling, or bronzewing barbules. See Section 4.2 of our paper
/// @param melaninRatio: a spatially varying quantity, modeling the percentage of melanin in regions where melanosomes are dispersed in keratin or air
void Geometry::simpleLayerStructures(MatrixXd melaninRatio) {
    parallel_for(nPoints, [&](int start, int end) {
    for (int i = start; i < end; i++) {
        Vector2d prev(xtop(i), ytop(i)), curr(xtop(i + 1), ytop(i + 1)), next(xtop(i + 2), ytop(i + 2));
        Vector2d t = (next - prev) / (next - prev).norm();
        double a = t(0), b = t(1);
        double c = -a * curr(0) - b * curr(1);
        vector<Vector2d> pts;
        for (int j = 1; j < nComponents; j++) {
            VectorXd xvec = xyvals.block(info(j, 0), 0, info(j, 1), 1);
            VectorXd yvec = xyvals.block(info(j, 0), 1, info(j, 1), 1);
            Vector2d pt = intersectLayer(a, b, c, xvec, yvec);
            pts.push_back(pt);
        }
        layer_height(0, i) = sqrt(pow(pts[0].x() - curr(0), 2) + pow(pts[0].y() - curr(1), 2));
        for (int j = 1; j < pts.size(); j++)
            layer_height(j, i) = sqrt(pow(pts[j].x() - pts[j - 1].x(), 2) + pow(pts[j].y() - pts[j - 1].y(), 2));
        double ratio = melaninRatio(i);
        if (type == 1) {
            layer_ior(0, i) = air;
            layer_ior(1, i) = keratin;
            layer_ior(2, i) = air;
            layer_ior(3, i) = (1 - ratio) * air + ratio * melanin;
            layer_ior(4, i) = keratin;
            layer_ior(5, i) = air;
        } else if (type == 2) {
            layer_ior(0, i) = air;
            layer_ior(1, i) = keratin;
            layer_ior(2, i) = melanin;
            layer_ior(3, i) = (1 - ratio) * keratin + ratio * melanin;
            layer_ior(4, i) = melanin;
            layer_ior(5, i) = keratin;
            layer_ior(6, i) = air;
        } else {
            for (int n = 0; n < nComponents / 2; n++) {
                if (n == 0) {
                    layer_ior(n, i) = air;
                    layer_ior(nComponents - n, i) = air;
                } else if (n % 2 == 1) {
                    layer_ior(n, i) = keratin;
                    layer_ior(nComponents - n, i) = keratin;
                } else {
                    layer_ior(n, i) = melanin;
                    layer_ior(nComponents - n, i) = melanin;
                }
            }
            layer_ior(nComponents / 2, i) = (1 - ratio) * keratin + ratio * melanin;
        }
    }
    } );
}

/// @brief Analyze the spatially varying multilayer structures in hummingbird barbules (with pancake-like structures). See Section 4.2 of our paper
void Geometry::hummingbirdStructures() {
    parallel_for(nPoints, [&](int start, int end) {
    for (int i = start; i < end; i++) {
        Vector2d prev(xtop(i), ytop(i)), curr(xtop(i + 1), ytop(i + 1)), next(xtop(i + 2), ytop(i + 2));
        Vector2d t = (next - prev) / (next - prev).norm();
        double a = t(0), b = t(1);
        double c = -a * curr(0) - b * curr(1);
        vector<Vector2d> pts;
        layer_ior(0, i) = air;
        layer_ior(1, i) = keratin;
        int count = 2;
        // Intersect with the pancake structures
        int nMelanin = (nComponents - 2) / 4;
        for (int j = 0; j < nMelanin; j++) {
            VectorXd x1 = xyvals.block(info(4 * j + 1, 0), 0, info(4 * j + 1, 1), 1);
            VectorXd y1 = xyvals.block(info(4 * j + 1, 0), 1, info(4 * j + 1, 1), 1);
            VectorXd x2 = xyvals.block(info(4 * j + 2, 0), 0, info(4 * j + 2, 1), 1);
            VectorXd y2 = xyvals.block(info(4 * j + 2, 0), 1, info(4 * j + 2, 1), 1);
            VectorXd x3 = xyvals.block(info(4 * j + 3, 0), 0, info(4 * j + 3, 1), 1);
            VectorXd y3 = xyvals.block(info(4 * j + 3, 0), 1, info(4 * j + 3, 1), 1);
            VectorXd x4 = xyvals.block(info(4 * j + 4, 0), 0, info(4 * j + 4, 1), 1);
            VectorXd y4 = xyvals.block(info(4 * j + 4, 0), 1, info(4 * j + 4, 1), 1);
            double v1 = a * x1(0) + b * y1(0) + c;
            double v2 = a * x1(x1.rows() - 1) + b * y1(y1.rows() - 1) + c;
            double v7 = a * x4(0) + b * y4(0) + c;
            double v8 = a * x4(x4.rows() - 1) + b * y4(y4.rows() - 1) + c;
            if (v1 * v2 >= 0 || v7 * v8 >= 0)
                continue;
            Vector2d pt1 = intersectLayer(a, b, c, x1, y1);
            Vector2d pt4 = intersectLayer(a, b, c, x4, y4);
            double v3 = a * x2(0) + b * y2(0) + c;
            double v4 = a * x2(x2.rows() - 1) + b * y2(y2.rows() - 1) + c;
            double v5 = a * x3(0) + b * y3(0) + c;
            double v6 = a * x3(x3.rows() - 1) + b * y3(y3.rows() - 1) + c;
            if (v3 * v4 >= 0 || v5 * v6 >= 0) {
                pts.push_back(pt1);
                pts.push_back(pt4);
                layer_ior(count, i) = melanin;
                layer_ior(count + 1, i) = keratin;
                count = count + 2;
            } else {
                Vector2d pt2 = intersectLayer(a, b, c, x2, y2);
                Vector2d pt3 = intersectLayer(a, b, c, x3, y3);
                pts.push_back(pt1);
                pts.push_back(pt2);
                pts.push_back(pt3);
                pts.push_back(pt4);
                layer_ior(count, i) = melanin;
                layer_ior(count + 1, i) = air;
                layer_ior(count + 2, i) = melanin;
                layer_ior(count + 3, i) = keratin;
                count = count + 4;
            }    
        }
        // Intersect with the bottom boundary
        VectorXd xvec = xyvals.block(info(nComponents - 1, 0), 0, info(nComponents - 1, 1), 1);
        VectorXd yvec = xyvals.block(info(nComponents - 1, 0), 1, info(nComponents - 1, 1), 1);
        Vector2d pt = intersectLayer(a, b, c, xvec, yvec);
        pts.push_back(pt);
        layer_ior(count, i) = air;
        layer_height(0, i) = sqrt(pow(pts[0].x() - curr(0), 2) + pow(pts[0].y() - curr(1), 2));
        for (int j = 1; j < pts.size(); j++)
            layer_height(j, i) = sqrt(pow(pts[j].x() - pts[j - 1].x(), 2) + pow(pts[j].y() - pts[j - 1].y(), 2));
    }
    } );
}

/// @brief Find the intersection point of a given line with one section (layer) in a barbule structure
/// @param a, b, c: specifying the line as ax + by + c = 0
/// @param xvec, yvec: a specific layer boundary in a barbule structure
Vector2d Geometry::intersectLayer(double a, double b, double c, VectorXd xvec, VectorXd yvec) {
    int num = xvec.rows();
    double v1 = a * xvec(0) + b * yvec(0) + c;
    double v2 = a * xvec(num - 1) + b * yvec(num - 1) + c;
    double x1, y1, x2, y2;
    if (v1 * v2 >= 0) {
        if (abs(v1) <= abs(v2)) {
            x1 = xvec(0);
            y1 = yvec(0);
            x2 = xvec(1);
            y2 = yvec(1);
        } else {
            x1 = xvec(num - 2);
            y1 = yvec(num - 2);
            x2 = xvec(num - 1);
            y2 = yvec(num - 1);
        }
    } else {
        int ind1 = 0, ind2 = num - 1;
        while (ind2 - ind1 > 1) {
            int ind = (ind1 + ind2) / 2;
            double v = a * xvec(ind) + b * yvec(ind) + c;
            if (abs(v) < 1e-4) {
                Vector2d pt(xvec(ind), yvec(ind));
                return pt;
            }
            if (v1 * v > 0) {
                ind1 = ind;
                v1 = v;
            } else {
                ind2 = ind;
                v2 = v;
            }
        }
        x1 = xvec(ind1);
        y1 = yvec(ind1);
        x2 = xvec(ind2);
        y2 = yvec(ind2);
    }

    // Compute the intersection point
    Vector2d pt1(x1, y1), pt2(x2, y2);
    Vector2d t = (pt2 - pt1) / (pt2 - pt1).norm();
    double a0 = -t(1), b0 = t(0), c0 = -a0 * x1 - b0 * y1;
    double x = (c0 * b - c * b0) / (a * b0 - a0 * b);
    double y = (c0 * a - c * a0) / (b * a0 - b0 * a);
    Vector2d pt(x, y);
    return pt;
}

/// @brief Fill in the photonic crystal reflection coefficients for mallard, magpie, or peacock barbules. See Section 4.3 of our paper
/// @param coef: the raw input reflection coefficients computed from FDTD simulations using the MEEP software
void Geometry::photonicReference(MatrixXd coef) {
    VectorXd cos_theta = VectorXd::LinSpaced(nTheta, 1.0 / nTheta, 1.0);
    VectorXd cos_phi = VectorXd::LinSpaced(nPhi, 0.0, 1.0);

    // Construct the coefficient table
    for (int i = 0; i < nTheta; i++) {
        vector<MatrixXcd> theta_coef;
        int n_bound = ceil(2.0 * cos_theta(i) * period / lambda);
        for (int n = 0; n <= n_bound; n++) {
            MatrixXcd n_coef = MatrixXcd::Zero(nPhi, 13);
            for (int j = 1; j < nPhi; j++) {
                n_coef(j, 0) = cos_phi(j);
                double sin_r = -sqrt(1 - cos_phi(j) * cos_phi(j)) + n * lambda / (cos_theta(i) * period);
                if (abs(sin_r) <= 1.0) {
                    int index = nTheta * (nPhi - 1) * n + (nPhi - 1) * i + (j - 1);
                    for (int col = 0; col < 12; col++)
                        n_coef(j, col + 1) = coef(col, 2 * index) + cunit * coef(col, 2 * index + 1);
                }
            }

            // Filter out missing values
            vector<int> valid_inds;
            for (int j = 0; j < nPhi; j++) {
                if (j == 0 || !isEmptySim(n_coef.block(j, 1, 1, 12)))
                    valid_inds.push_back(j);
            }
            if (valid_inds.size() <= 1)
                break;
            MatrixXcd phi_coef = MatrixXcd::Zero(valid_inds.size(), 13);
            for (int j = 0; j < valid_inds.size(); j++)
                phi_coef.block(j, 0, 1, 13) = n_coef.block(valid_inds[j], 0, 1, 13);

            // Add in coefficients for grazing incidence (phi = 90 degree)
            if (n == 0) {
                phi_coef(0, 1) = sqrt(1 - cos_theta(i) * cos_theta(i));
                phi_coef(0, 3) = -cos_theta(i);
                phi_coef(0, 5) = -1.0;
                phi_coef(0, 8) = 1.0;
                phi_coef(0, 10) = sqrt(1 - cos_theta(i) * cos_theta(i));
                phi_coef(0, 12) = -cos_theta(i);
            }
            theta_coef.push_back(phi_coef);
        }
        photonic_coef.push_back(theta_coef);
    }
}

/// @brief Helper function that tests if a given row of data is all zero-valued
/// @param row_data: the given row of data
bool Geometry::isEmptySim(MatrixXcd row_data) {
    for (int check = 0; check < 12; check++) {
        if (real(row_data(0, check)) != 0 || imag(row_data(0, check)) != 0)
            return false;
    }
    return true;
}