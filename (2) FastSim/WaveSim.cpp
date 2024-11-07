/* This module implements our approximate (but very accurate) wave simulations on different types of barbules.
Users are welcome to read but do not need to fully understand methods in this file.
PLEASE DO NOT MODIFY THE FOLLOWING CODE. */

#include "WaveSim.h"

/// @brief Some important initializations
/// @param geom: pointer to a Geometry object describing a barbule geometry
WaveSim::WaveSim(Geometry* geom) {
    this->geom = geom;
    type = geom->type;
    lambda = geom->lambda;
    air = 1.0;
    k = 2.0 * M_PI * air / lambda;
    omega = 2.0 * M_PI * c0 / lambda;
    eps = 1 / (mu * c0 * c0) * air * air;
    if (type >= 5)
        nTheta = geom->nTheta;

    // Gaussian quadrature for numerical integration
    order = 3;
    pvals = quadrature_points.block(order - 1, 0, 1, order).transpose();
    wvals = quadrature_weights.block(order - 1, 0, 1, order).transpose();
}

/// @brief The umbrella function that sets the simulation running, for computing the incoming light's incident power and the scattering quantities
/// @param theta, phi: specifying the direction of the incident light. See Section 4.4 of our paper
/// @param resolution1, resolution2: specified resolutions used for computing the scattering quantities
/// @param regularity: a spatially varying quantity, modeling the level of organization and regularity in mallard, magpie, or peacock barbules
void WaveSim::simulate(double theta, double phi, int resolution1, int resolution2, MatrixXd regularity) {
    xtop = geom->xtop;
    ytop = geom->ytop;
    computeIncidentPower(theta, phi);
    computeSurfaceCurrents(theta, phi, regularity);
    computeFarField(theta, resolution1, false, farfield1);
    computeFarField(theta, resolution2, true, farfield2);
}

/// @brief Compute the incident field's total incident power on the barbule, by integrating the irradiance over the barbule top surface
/// @param theta, phi: specifying the direction of the incident light. See Section 4.4 of our paper
void WaveSim::computeIncidentPower(double theta, double phi) {
    double kx = -k * cos(theta) * sin(phi);
    double ky = -k * cos(theta) * cos(phi);
    inci_p1 = 0;
    inci_p2 = 0;
    double range = xtop.maxCoeff() - xtop.minCoeff();
    double sd = 0.25 * range;

    // Integrating the irradiance on the barbule top surface
    for (int i = 0; i < xtop.rows() - 2; i++) {
        Vector2d prev(xtop(i), ytop(i)), curr(xtop(i + 1), ytop(i + 1)), next(xtop(i + 2), ytop(i + 2));
        double scale = exp(-curr(0) * curr(0) / (sd * sd));
        Vector2d t = (next - prev) / (next - prev).norm();
        Vector3d n(-t(1), t(0), 0);
        double d = 0.5 * (curr - prev).norm() + 0.5 * (next - curr).norm();
        dcomp phase = exp(cunit * (kx * curr(0) + ky * curr(1)));
        Vector3d poynting1, poynting2;
        dcomp Ex1 = -sin(theta) * sin(phi) * phase;
        dcomp Ey1 = -sin(theta) * cos(phi) * phase;
        dcomp Ez1 = cos(theta) * phase;
        dcomp Hx1 = -air / eta0 * cos(phi) * phase;
        dcomp Hy1 = air / eta0 * sin(phi) * phase;
        dcomp Hz1 = 0;
        poynting1(0) = 0.5 * (Ey1 * conj(Hz1) - Ez1 * conj(Hy1)).real();
        poynting1(1) = 0.5 * (Ez1 * conj(Hx1) - Ex1 * conj(Hz1)).real();
        poynting1(2) = 0.5 * (Ex1 * conj(Hy1) - Ey1 * conj(Hx1)).real();
        dcomp Ex2 = cos(phi) * phase;
        dcomp Ey2 = -sin(phi) * phase;
        dcomp Ez2 = 0;
        dcomp Hx2 = -air / eta0 * sin(theta) * sin(phi) * phase;
        dcomp Hy2 = -air / eta0 * sin(theta) * cos(phi) * phase;
        dcomp Hz2 = air / eta0 * cos(theta) * phase;
        poynting2(0) = 0.5 * (Ey2 * conj(Hz2) - Ez2 * conj(Hy2)).real();
        poynting2(1) = 0.5 * (Ez2 * conj(Hx2) - Ex2 * conj(Hz2)).real();
        poynting2(2) = 0.5 * (Ex2 * conj(Hy2) - Ey2 * conj(Hx2)).real();
        if (n.dot(poynting1) < 0)
            inci_p1 += (float)(d * scale * scale * abs(n.dot(poynting1)));
        if (n.dot(poynting2) < 0)
            inci_p2 += (float)(d * scale * scale * abs(n.dot(poynting2)));
    }
}

/// @brief Approximate the fictitious current densities on the barbule top surface, which will be used for computing far field scattering
/// @param theta, phi: specifying the direction of the incident light. See Section 4.4 of our paper
/// @param regularity: a spatially varying quantity, modeling the level of organization and regularity in mallard, magpie, or peacock barbules
void WaveSim::computeSurfaceCurrents(double theta, double phi, MatrixXd regularity) {
    Vector3d ei(cos(theta) * sin(phi), cos(theta) * cos(phi), sin(theta));
    Vector3d ki = -k * ei, kt = -k * ei;
    double range = xtop.maxCoeff() - xtop.minCoeff();
    double sd = 0.25 * range;

    // Estimate current densities at each point on the barbule
    Jt1 = VectorXcf::Zero(xtop.rows() - 2);
    Jz1 = VectorXcf::Zero(xtop.rows() - 2);
    Mt1 = VectorXcf::Zero(xtop.rows() - 2);
    Mz1 = VectorXcf::Zero(xtop.rows() - 2);
    Jt2 = VectorXcf::Zero(xtop.rows() - 2);
    Jz2 = VectorXcf::Zero(xtop.rows() - 2);
    Mt2 = VectorXcf::Zero(xtop.rows() - 2);
    Mz2 = VectorXcf::Zero(xtop.rows() - 2);
    parallel_for(xtop.rows() - 2, [&](int start, int end) {
    for (int i = start; i < end; i++) {
        Vector2d prev(xtop(i), ytop(i)), curr(xtop(i + 1), ytop(i + 1)), next(xtop(i + 2), ytop(i + 2));
        float scale = exp(-curr(0) * curr(0) / (sd * sd));
        if (type >= 5)
            scale = scale * regularity(i);
        Vector2d t = (next - prev) / (next - prev).norm();
        Vector3d n(-t(1), t(0), 0);
        dcomp phase = exp(cunit * (ki(0) * curr(0) + ki(1) * curr(1)));
        double tilt = acos(ei(0) * n(0) + ei(1) * n(1) + ei(2) * n(2));
        if (tilt >= M_PI / 2 - 1e-3)
            continue;
        Vector3d er = 2 * cos(tilt) * n - ei;
        Vector3d kr = k * er;
        Vector3d vPerp = ei.cross(n);
        if (vPerp.norm() < 1e-4)
            vPerp << 0, 0, 1;
        else
            vPerp = vPerp / vPerp.norm();

        // Incident field from two polarizations
        double Ex1 = -sin(theta) * sin(phi);
        double Ey1 = -sin(theta) * cos(phi);
        double Ez1 = cos(theta);
        double Hx1 = -air / eta0 * cos(phi);
        double Hy1 = air / eta0 * sin(phi);
        double Hz1 = 0;
        double c_TE1 = Ex1 * vPerp(0) + Ey1 * vPerp(1) + Ez1 * vPerp(2);
        double c_TM1 = Hx1 * vPerp(0) + Hy1 * vPerp(1) + Hz1 * vPerp(2);
        double Ex2 = cos(phi);
        double Ey2 = -sin(phi);
        double Ez2 = 0;
        double Hx2 = -air / eta0 * sin(theta) * sin(phi);
        double Hy2 = -air / eta0 * sin(theta) * cos(phi);
        double Hz2 = air / eta0 * cos(theta);
        double c_TE2 = Ex2 * vPerp(0) + Ey2 * vPerp(1) + Ez2 * vPerp(2);
        double c_TM2 = Hx2 * vPerp(0) + Hy2 * vPerp(1) + Hz2 * vPerp(2);

        // Compute reflected field and current density
        Vector3cd Er1(0, 0, 0), Hr1(0, 0, 0), Er2(0, 0, 0), Hr2(0, 0, 0);
        if (type >= 5)
            queryPhotonics(theta, phi, curr, t, Er1, Hr1, Er2, Hr2);
        else {
            int nH = 0;
            while (nH < geom->layer_height.rows() && geom->layer_height(nH, i) != 0)
                nH = nH + 1;
            VectorXd height = geom->layer_height.block(0, i, nH, 1);
            VectorXcd ior = geom->layer_ior.block(0, i, nH + 2, 1);
            Vector2cd gamma = compositeR(height, ior, tilt);
            Er1 = gamma(0) * c_TE1 * phase * vPerp - gamma(1) * c_TM1 * phase * kr.cross(vPerp) / (omega * eps);
            Hr1 = gamma(0) * c_TE1 * phase * kr.cross(vPerp) / (omega * mu) + gamma(1) * c_TM1 * phase * vPerp;
            Er2 = gamma(0) * c_TE2 * phase * vPerp - gamma(1) * c_TM2 * phase * kr.cross(vPerp) / (omega * eps);
            Hr2 = gamma(0) * c_TE2 * phase * kr.cross(vPerp) / (omega * mu) + gamma(1) * c_TM2 * phase * vPerp;
        }
        Jt1(i) = scale * (fcomp)Hr1(2);
        Jz1(i) = scale * (fcomp)(n(0) * Hr1(1) - n(1) * Hr1(0));
        Mt1(i) = scale * (fcomp)(-Er1(2) / eta0);
        Mz1(i) = scale * (fcomp)((n(1) * Er1(0) - n(0) * Er1(1)) / eta0);
        Jt2(i) = scale * (fcomp)Hr2(2);
        Jz2(i) = scale * (fcomp)(n(0) * Hr2(1) - n(1) * Hr2(0));
        Mt2(i) = scale * (fcomp)(-Er2(2) / eta0);
        Mz2(i) = scale * (fcomp)((n(1) * Er2(0) - n(0) * Er2(1)) / eta0);
    }
    } );
}

/// @brief Compute the composite reflectivity of an ideal multilayer structure
/// @param height: a vector containing the thicknesses of each layer in the ideal multilayer structure
/// @param ior: a vector containing the IORs of each layer in the ideal multilayer structure
/// @param tilt: the incident angle
Vector2cd WaveSim::compositeR(VectorXd height, VectorXcd ior, double tilt) {
    int nBoundary = ior.rows() - 1;
    VectorXcd cos_theta = VectorXcd::Zero(nBoundary + 1);
    for (int i = 0; i < nBoundary + 1; i++) {
        dcomp sin_medium = sin(tilt) / ior(i);
        cos_theta(i) = sqrt(1.0 - sin_medium * sin_medium);
    }

    // For TE polarization, we compute r and t with respect to Ez; for TM polarization, we compute with respect to Hz
    dcomp gamma_te = (ior(nBoundary - 1) * cos_theta(nBoundary - 1) - ior(nBoundary) * cos_theta(nBoundary)) / (ior(nBoundary - 1) * cos_theta(nBoundary - 1) + ior(nBoundary) * cos_theta(nBoundary));
    dcomp gamma_tm = (ior(nBoundary) * cos_theta(nBoundary - 1) - ior(nBoundary - 1) * cos_theta(nBoundary)) / (ior(nBoundary) * cos_theta(nBoundary - 1) + ior(nBoundary - 1) * cos_theta(nBoundary));
    for (int i = nBoundary - 1; i >= 1; i--) {
        dcomp r_te = (ior(i - 1) * cos_theta(i - 1) - ior(i) * cos_theta(i)) / (ior(i - 1) * cos_theta(i - 1) + ior(i) * cos_theta(i));
        dcomp r_tm = (ior(i) * cos_theta(i - 1) - ior(i - 1) * cos_theta(i)) / (ior(i) * cos_theta(i - 1) + ior(i - 1) * cos_theta(i));
        dcomp phase = exp(cunit * k * ior(i) * height(i - 1) * cos_theta(i));
        gamma_te = (r_te + gamma_te * phase * phase) / (1.0 + r_te * gamma_te * phase * phase);
        gamma_tm = (r_tm + gamma_tm * phase * phase) / (1.0 + r_tm * gamma_tm * phase * phase);
    }
    Vector2cd gamma(gamma_te, gamma_tm);
    return gamma;
}

/// @brief Estimate the reflected field value at a given point along the surface of a photonic crystal type barbule
/// @param theta, phi: specifying the direction of the incident light. See Section 4.4 of our paper
/// @param pt: the point at which we want to estimate the reflected field values
/// @param tvec: the local tangent direction of the barbule top surface
/// @param Er1, Hr1, Er2, Hr2: the approximated reflected eletric and magnetic fields, for two different polarizations
void WaveSim::queryPhotonics(double theta, double phi, Vector2d pt, Vector2d tvec, Vector3cd& Er1, Vector3cd& Hr1, Vector3cd& Er2, Vector3cd& Hr2) {
    double tilt = asin(tvec(1));
    double cos_theta = cos(theta), cos_phi = cos(phi + tilt);
    if (cos_phi <= 0)
        return;
    vector<MatrixXcd> coef_repo = geom->photonic_coef[round(cos_theta * nTheta) - 1];

    // Compute the reflected field components and coefficients at the current surface point
    vector<double> phi_r;
    vector<VectorXcd> coefs;
    for (int n = 0; n < coef_repo.size(); n++) {
        double sin_r = -sin(phi + tilt) + n * lambda / (cos_theta * geom->period);
        if (sin(phi + tilt) < 0)
            sin_r = -sin(phi + tilt) - n * lambda / (cos_theta * geom->period);
        if (abs(sin_r) <= 1) {
            phi_r.push_back(asin(sin_r));
            MatrixXcd n_coef = coef_repo[n];
            int ind = 0;
            while (ind < n_coef.rows() && real(n_coef(ind, 0)) <= cos_phi)
                ind = ind + 1;
            if (ind == 0)
                ind = 1;
            if (ind >= n_coef.rows())
                ind = n_coef.rows() - 1;
            double r = (cos_phi - real(n_coef(ind - 1, 0))) / (real(n_coef(ind, 0)) - real(n_coef(ind - 1, 0)));
            MatrixXcd interp = (1.0 - r) * n_coef.block(ind - 1, 1, 1, 12) + r * n_coef.block(ind, 1, 1, 12);
            VectorXcd coefEH = interp.transpose();
            if (sin(phi + tilt) < 0) {
                coefEH(0) = -coefEH(0);
                coefEH(4) = -coefEH(4);
                coefEH(5) = -coefEH(5);
                coefEH(7) = -coefEH(7);
                coefEH(8) = -coefEH(8);
                coefEH(9) = -coefEH(9);
            }
            coefs.push_back(coefEH);
        }
    }

    // Compute the reflected field values at the current surface point
    double kx = -k * cos(theta) * sin(phi);
    double ky = -k * cos(theta) * cos(phi);
    dcomp phase = exp(cunit * (kx * pt(0) + ky * pt(1)));
    for (int n = 0; n < phi_r.size(); n++) {
        VectorXcd coefEH = coefs[n];
        Er1(0) = Er1(0) + (cos(tilt) * coefEH(0) - sin(tilt) * coefEH(1)) * phase;
        Er1(1) = Er1(1) + (sin(tilt) * coefEH(0) + cos(tilt) * coefEH(1)) * phase;
        Er1(2) = Er1(2) + coefEH(2) * phase;
        Hr1(0) = Hr1(0) + air / eta0 * (cos(tilt) * coefEH(3) - sin(tilt) * coefEH(4)) * phase;
        Hr1(1) = Hr1(1) + air / eta0 * (sin(tilt) * coefEH(3) + cos(tilt) * coefEH(4)) * phase;
        Hr1(2) = Hr1(2) + air / eta0 * coefEH(5) * phase;
        Er2(0) = Er2(0) + (cos(tilt) * coefEH(6) - sin(tilt) * coefEH(7)) * phase;
        Er2(1) = Er2(1) + (sin(tilt) * coefEH(6) + cos(tilt) * coefEH(7)) * phase;
        Er2(2) = Er2(2) + coefEH(8) * phase;
        Hr2(0) = Hr2(0) + air / eta0 * (cos(tilt) * coefEH(9) - sin(tilt) * coefEH(10)) * phase;
        Hr2(1) = Hr2(1) + air / eta0 * (sin(tilt) * coefEH(9) + cos(tilt) * coefEH(10)) * phase;
        Hr2(2) = Hr2(2) + air / eta0 * coefEH(11) * phase;
    }
}

/// @brief Compute the far field scattering quantities
/// @param theta: specifying the longitudinal angle of the incident light. See Section 4.4 of our paper
/// @param resolution: specified resolution used for computing the scattering quantities
/// @param use_projected: helper argument that determines how the outgoing azimuthal angles are selected
/// @param farfield: the far field scattering quantities
void WaveSim::computeFarField(double theta, int resolution, bool use_projected, MatrixXf& farfield) {
    farfield = MatrixXf::Zero(resolution, 14);
    parallel_for(resolution, [&](int start, int end) {
    for (int count = start; count < end; count++) {
        double phi = -0.5 * M_PI + 0.5 * M_PI / resolution + count * M_PI / resolution;
        if (use_projected) {
            double x0 = -1.0 + 1.0 / resolution + 2.0 * count / resolution;
            if (abs(x0) > cos(theta))
                continue;
            phi = asin(x0 / cos(theta));
        }

        // Polarization 1
        MatrixXcf field1 = evaluateScattering(xtop, ytop, theta, phi, Jt1, Jz1, Mt1, Mz1);
        Vector3f p1;
        p1(0) = 0.5f * real(field1(1, 0) * conj(field1(2, 1)) - field1(2, 0) * conj(field1(1, 1)));
        p1(1) = 0.5f * real(field1(2, 0) * conj(field1(0, 1)) - field1(0, 0) * conj(field1(2, 1)));
        p1(2) = 0.5f * real(field1(0, 0) * conj(field1(1, 1)) - field1(1, 0) * conj(field1(0, 1)));
        farfield(count, 0) = p1.norm() / inci_p1;
        farfield(count, 1) = field1(0, 0).real() / sqrt(inci_p1);
        farfield(count, 2) = field1(0, 0).imag() / sqrt(inci_p1);
        farfield(count, 3) = field1(1, 0).real() / sqrt(inci_p1);
        farfield(count, 4) = field1(1, 0).imag() / sqrt(inci_p1);
        farfield(count, 5) = field1(2, 0).real() / sqrt(inci_p1);
        farfield(count, 6) = field1(2, 0).imag() / sqrt(inci_p1);

        // Polarization 2
        MatrixXcf field2 = evaluateScattering(xtop, ytop, theta, phi, Jt2, Jz2, Mt2, Mz2);
        Vector3f p2;
        p2(0) = 0.5f * real(field2(1, 0) * conj(field2(2, 1)) - field2(2, 0) * conj(field2(1, 1)));
        p2(1) = 0.5f * real(field2(2, 0) * conj(field2(0, 1)) - field2(0, 0) * conj(field2(2, 1)));
        p2(2) = 0.5f * real(field2(0, 0) * conj(field2(1, 1)) - field2(1, 0) * conj(field2(0, 1)));
        farfield(count, 7) = p2.norm() / inci_p2;
        farfield(count, 8) = (float)eta0 * field2(0, 1).real() / sqrt(inci_p2);
        farfield(count, 9) = (float)eta0 * field2(0, 1).imag() / sqrt(inci_p2);
        farfield(count, 10) = (float)eta0 * field2(1, 1).real() / sqrt(inci_p2);
        farfield(count, 11) = (float)eta0 * field2(1, 1).imag() / sqrt(inci_p2);
        farfield(count, 12) = (float)eta0 * field2(2, 1).real() / sqrt(inci_p2);
        farfield(count, 13) = (float)eta0 * field2(2, 1).imag() / sqrt(inci_p2);
    }
    } );
}

/// @brief Evaluate the scattered field quantities along a specified outgoing direction
/// @param xvals, yvals: the 2D geometry of the barbule top surface
/// @param theta, phi: specifying the queried outgoing direction
/// @param Jt, Jz, Mt, Mz: the fictitious current densities along the barbule top surface
MatrixXcf WaveSim::evaluateScattering(VectorXd xvals, VectorXd yvals, double theta, double phi, VectorXcf Jt, VectorXcf Jz, VectorXcf Mt, VectorXcf Mz) {
    Vector3d rho(sin(phi), cos(phi), 0), z(0, 0, 1), normal(cos(phi), -sin(phi), 0);
    Vector3cf Es(0, 0, 0), Hs(0, 0, 0);
    int npts = xvals.rows() - 2;
    
    // Accumulate the contribution of each basis element into the far field quantities
    for (int i = 0; i < npts; i++) {
        Vector3d prev(xvals(i), yvals(i), 0), curr(xvals(i + 1), yvals(i + 1), 0), next(xvals(i + 2), yvals(i + 2), 0);
        double delta1 = (curr - prev).norm();
        double delta2 = (next - curr).norm();
        Vector3d t1 = (curr - prev) / delta1;
        Vector3d t2 = (next - curr) / delta2;
        Vector3d n1(-t1(1), t1(0), 0);
        Vector3d n2(-t2(1), t2(0), 0);
        fcomp jt = Jt(i), jz = Jz(i), mt = Mt(i) * (float)eta0, mz = Mz(i) * (float)eta0;
        
        // First half of the basis function
        for (int count = 0; count < order; count++) {
            double r = 0.5 * pvals(count) + 0.5;
            Vector3d pt = (1 - r) * prev + r * curr;
            dcomp g0 = sqrt(2.0 / (cunit * M_PI * k * cos(theta))) * exp(-cunit * k * cos(theta) * (sin(phi) * pt(0) + cos(phi) * pt(1)));
            dcomp g1 = -g0 * cunit;
            Vector3cd field1 = wvals(count) * delta1 / 2.0 * g0 * r * t1;
            Vector3cd field2 = wvals(count) / 2.0 * g1 * rho;
            Vector3cd field3 = wvals(count) / 2.0 * g0 * z;
            Vector3cd field4 = wvals(count) * delta1 / 2.0 * g0 * r * z;
            Vector3cd field5 = wvals(count) * delta1 / 2.0 * g1 * r * rho;
            Vector3cd field6 = wvals(count) * delta1 / 2.0 * g1 * r * (rho(0) * t1(1) - rho(1) * t1(0)) * z;
            Vector3cd field7 = wvals(count) * delta1 / 2.0 * g0 * r * n1;
            Vector3cd field8 = wvals(count) * delta1 / 2.0 * g1 * r * normal;
            Vector3cd Es1 = -omega * mu / 4.0 * field1 + k * cos(theta) / (4.0 * omega * eps) * field2 + cunit * k * sin(theta) / (4.0 * omega * eps) * field3;
            Vector3cd Es2 = (-omega * mu / 4.0 + k * k * sin(theta) * sin(theta) / (4.0 * omega * eps)) * field4 - cunit * k * k * cos(theta) * sin(theta) / (4.0 * omega * eps) * field5;
            Vector3cd Es3 = cunit * k * cos(theta) / 4.0 * field6 - k * sin(theta) / 4.0 * field7;
            Vector3cd Es4 = cunit * k * cos(theta) / 4.0 * field8;
            Es += jt * Es1.cast<fcomp>() + jz * Es2.cast<fcomp>() + mt * Es3.cast<fcomp>() + mz * Es4.cast<fcomp>();
            Vector3cd Hs1 = -cunit * k * cos(theta) / 4.0 * field6 + k * sin(theta) / 4.0 * field7;
            Vector3cd Hs2 = -cunit * k * cos(theta) / 4.0 * field8;
            Vector3cd Hs3 = -omega * eps / 4.0 * field1 + k * cos(theta) / (4.0 * omega * mu) * field2 + cunit * k * sin(theta) / (4.0 * omega * mu) * field3;
            Vector3cd Hs4 = (-omega * eps / 4.0 + k * k * sin(theta) * sin(theta) / (4.0 * omega * mu)) * field4 - cunit * k * k * cos(theta) * sin(theta) / (4.0 * omega * mu) * field5;
            Hs += jt * Hs1.cast<fcomp>() + jz * Hs2.cast<fcomp>() + mt * Hs3.cast<fcomp>() + mz * Hs4.cast<fcomp>();
        }
        
        // Second half of the basis function
        for (int count = 0; count < order; count++) {
            double r = 0.5 * pvals(count) + 0.5;
            Vector3d pt = (1 - r) * curr + r * next;
            dcomp g0 = sqrt(2.0 / (cunit * M_PI * k * cos(theta))) * exp(-cunit * k * cos(theta) * (sin(phi) * pt(0) + cos(phi) * pt(1)));
            dcomp g1 = -g0 * cunit;
            Vector3cd field1 = wvals(count) * delta2 / 2.0 * g0 * (1 - r) * t2;
            Vector3cd field2 = -wvals(count) / 2.0 * g1 * rho;
            Vector3cd field3 = -wvals(count) / 2.0 * g0 * z;
            Vector3cd field4 = wvals(count) * delta2 / 2.0 * g0 * (1 - r) * z;
            Vector3cd field5 = wvals(count) * delta2 / 2.0 * g1 * (1 - r) * rho;
            Vector3cd field6 = wvals(count) * delta2 / 2.0 * g1 * (1 - r) * (rho(0) * t2(1) - rho(1) * t2(0)) * z;
            Vector3cd field7 = wvals(count) * delta2 / 2.0 * g0 * (1 - r) * n2;
            Vector3cd field8 = wvals(count) * delta2 / 2.0 * g1 * (1 - r) * normal;
            Vector3cd Es1 = -omega * mu / 4.0 * field1 + k * cos(theta) / (4.0 * omega * eps) * field2 + cunit * k * sin(theta) / (4.0 * omega * eps) * field3;
            Vector3cd Es2 = (-omega * mu / 4.0 + k * k * sin(theta) * sin(theta) / (4.0 * omega * eps)) * field4 - cunit * k * k * cos(theta) * sin(theta) / (4.0 * omega * eps) * field5;
            Vector3cd Es3 = cunit * k * cos(theta) / 4.0 * field6 - k * sin(theta) / 4.0 * field7;
            Vector3cd Es4 = cunit * k * cos(theta) / 4.0 * field8;
            Es += jt * Es1.cast<fcomp>() + jz * Es2.cast<fcomp>() + mt * Es3.cast<fcomp>() + mz * Es4.cast<fcomp>();
            Vector3cd Hs1 = -cunit * k * cos(theta) / 4.0 * field6 + k * sin(theta) / 4.0 * field7;
            Vector3cd Hs2 = -cunit * k * cos(theta) / 4.0 * field8;
            Vector3cd Hs3 = -omega * eps / 4.0 * field1 + k * cos(theta) / (4.0 * omega * mu) * field2 + cunit * k * sin(theta) / (4.0 * omega * mu) * field3;
            Vector3cd Hs4 = (-omega * eps / 4.0 + k * k * sin(theta) * sin(theta) / (4.0 * omega * mu)) * field4 - cunit * k * k * cos(theta) * sin(theta) / (4.0 * omega * mu) * field5;
            Hs += jt * Hs1.cast<fcomp>() + jz * Hs2.cast<fcomp>() + mt * Hs3.cast<fcomp>() + mz * Hs4.cast<fcomp>();
        }
    }
    MatrixXcf fields(3, 2);
    fields.block(0, 0, 3, 1) = Es;
    fields.block(0, 1, 3, 1) = Hs;
    return fields;
}