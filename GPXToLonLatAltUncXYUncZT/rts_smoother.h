#ifndef RTS_SMOOTHER_H
#define RTS_SMOOTHER_H
#include "Common.hpp"

struct GPSMeasurement {
    double lat;       // degrees
    double lon;       // degrees
    double alt;       // meters
    double horiz_std = 0.0; // horizontal position std dev (m)
    double vert_std  = 0.0; // vertical position std dev (m)
    double timestamp; // seconds
    bool   has_std = false;
};

struct State {
    std::array<double, 3> pos{}; // meters in ECEF
    std::array<double, 3> vel{}; // m/s in ECEF
};

using Matrix6 = std::array<std::array<double, 6>, 6>;
using Matrix3 = std::array<std::array<double, 3>, 3>;

struct StepData {
    State predicted_state;
    Matrix6 predicted_cov;
    State filtered_state;
    Matrix6 filtered_cov;
    double dt = 0.0;
	double timestamp = 0.0;
};

// Convert degrees to radians
 double deg2rad(double deg) { return deg * std::numbers::pi / 180.0; }
 double rad2deg(double rad) { return rad * 180.0 / std::numbers::pi; }

// Convert geodetic coordinates to ECEF (WGS84)
 std::array<double,3> geodeticToECEF(double lat_deg, double lon_deg, double alt_m) {
    const double a = 6378137.0; // WGS84 major axis
    const double e2 = 6.69437999014e-3;
    double lat = deg2rad(lat_deg);
    double lon = deg2rad(lon_deg);
    double N = a / std::sqrt(1.0 - e2 * std::sin(lat) * std::sin(lat));
    double x = (N + alt_m) * std::cos(lat) * std::cos(lon);
    double y = (N + alt_m) * std::cos(lat) * std::sin(lon);
    double z = ((1 - e2) * N + alt_m) * std::sin(lat);
    return {x, y, z};
}

// Convert ECEF coordinates to geodetic (WGS84)
 std::array<double,3> ecefToGeodetic(double x, double y, double z) {
    const double a = 6378137.0; // WGS84 major axis
    const double e2 = 6.69437999014e-3;
    const double b = a * std::sqrt(1.0 - e2);
    const double ep = std::sqrt((a*a - b*b) / (b*b));

    double p = std::sqrt(x*x + y*y);
    double th = std::atan2(a*z, b*p);
    double lon = std::atan2(y, x);
    double lat = std::atan2(z + ep*ep*b*std::pow(std::sin(th),3),
                            p - e2*a*std::pow(std::cos(th),3));
    double N = a / std::sqrt(1.0 - e2*std::sin(lat)*std::sin(lat));
    double alt = p / std::cos(lat) - N;

    return {rad2deg(lat), rad2deg(lon), alt};
}

 void estimateUncertainties(std::vector<GPSMeasurement>& meas) {
    if(meas.empty()) return;
    for(size_t i=0;i<meas.size();++i) {
        if(meas[i].has_std) continue;
        if(i==0 || i==meas.size()-1) {
            meas[i].horiz_std = 5.0;
            meas[i].vert_std = 5.0;
            meas[i].has_std = true;
            continue;
        }
        auto p_prev = geodeticToECEF(meas[i-1].lat, meas[i-1].lon, meas[i-1].alt);
        auto p_cur  = geodeticToECEF(meas[i].lat,  meas[i].lon,  meas[i].alt);
        auto p_next = geodeticToECEF(meas[i+1].lat, meas[i+1].lon, meas[i+1].alt);
        double dt_total = meas[i+1].timestamp - meas[i-1].timestamp;
        double dt_local = meas[i].timestamp - meas[i-1].timestamp;
        std::array<double,3> pred{};
        for(int k=0;k<3;++k) {
            double vel = (p_next[k]-p_prev[k])/dt_total;
            pred[k] = p_prev[k] + vel*dt_local;
        }
        std::array<double,3> diff{};
        for(int k=0;k<3;++k) diff[k] = p_cur[k] - pred[k];
        meas[i].horiz_std = std::sqrt(diff[0]*diff[0] + diff[1]*diff[1]);
        meas[i].vert_std = std::abs(diff[2]);
        meas[i].has_std = true;
    }
}

 Matrix3 inverse3(const Matrix3& m) {
    double a = m[0][0], b = m[0][1], c = m[0][2];
    double d = m[1][0], e = m[1][1], f = m[1][2];
    double g = m[2][0], h = m[2][1], i = m[2][2];
    double det = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    if(std::abs(det) < 1e-12) throw std::runtime_error("Matrix singular");
    Matrix3 inv;
    inv[0][0] =  (e*i - f*h)/det;
    inv[0][1] = -(b*i - c*h)/det;
    inv[0][2] =  (b*f - c*e)/det;
    inv[1][0] = -(d*i - f*g)/det;
    inv[1][1] =  (a*i - c*g)/det;
    inv[1][2] = -(a*f - c*d)/det;
    inv[2][0] =  (d*h - e*g)/det;
    inv[2][1] = -(a*h - b*g)/det;
    inv[2][2] =  (a*e - b*d)/det;
    return inv;
}

 Matrix6 identity6() {
    Matrix6 I{};
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) I[i][j] = (i==j) ? 1.0 : 0.0;
    return I;
}

 Matrix6 transpose6(const Matrix6& m) {
    Matrix6 t{};
    for(int i=0;i<6;++i)
        for(int j=0;j<6;++j)
            t[i][j] = m[j][i];
    return t;
}

 Matrix3 transpose3(const Matrix3& m) {
    Matrix3 t{};
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            t[i][j] = m[j][i];
    return t;
}

 Matrix6 multiply6(const Matrix6& a, const Matrix6& b) {
    Matrix6 r{};
    for(int i=0;i<6;++i)
        for(int j=0;j<6;++j) {
            double sum=0.0;
            for(int k=0;k<6;++k) sum += a[i][k]*b[k][j];
            r[i][j]=sum;
        }
    return r;
}

 std::array<double,6> multiply6(const Matrix6& a, const std::array<double,6>& v) {
    std::array<double,6> r{};
    for(int i=0;i<6;++i) {
        double sum=0.0;
        for(int k=0;k<6;++k) sum += a[i][k]*v[k];
        r[i]=sum;
    }
    return r;
}

 Matrix3 multiply3(const Matrix3& a, const Matrix3& b) {
    Matrix3 r{};
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j) {
            double sum=0.0;
            for(int k=0;k<3;++k) sum += a[i][k]*b[k][j];
            r[i][j]=sum;
        }
    return r;
}

 std::array<double,3> multiply3(const Matrix3& a, const std::array<double,3>& v) {
    std::array<double,3> r{};
    for(int i=0;i<3;++i) {
        double sum=0.0;
        for(int k=0;k<3;++k) sum += a[i][k]*v[k];
        r[i]=sum;
    }
    return r;
}

 Matrix6 add6(const Matrix6& a, const Matrix6& b) {
    Matrix6 r{};
    for(int i=0;i<6;++i)
        for(int j=0;j<6;++j)
            r[i][j] = a[i][j] + b[i][j];
    return r;
}

 Matrix6 subtract6(const Matrix6& a, const Matrix6& b) {
    Matrix6 r{};
    for(int i=0;i<6;++i)
        for(int j=0;j<6;++j)
            r[i][j] = a[i][j] - b[i][j];
    return r;
}

 Matrix3 add3(const Matrix3& a, const Matrix3& b) {
    Matrix3 r{};
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            r[i][j] = a[i][j] + b[i][j];
    return r;
}

 Matrix3 subtract3(const Matrix3& a, const Matrix3& b) {
    Matrix3 r{};
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            r[i][j] = a[i][j] - b[i][j];
    return r;
}

 std::array<double,6> subtract6(const std::array<double,6>& a, const std::array<double,6>& b) {
    std::array<double,6> r{};
    for(int i=0;i<6;++i) r[i]=a[i]-b[i];
    return r;
}

 std::array<double,6> add6(const std::array<double,6>& a, const std::array<double,6>& b) {
    std::array<double,6> r{};
    for(int i=0;i<6;++i) r[i]=a[i]+b[i];
    return r;
}


// Setup process noise Q for acceleration noise variance accel_var
 Matrix6 processNoise(double dt, double accel_var) {
    Matrix6 Q{};
    double dt2 = dt*dt;
    double dt3 = dt2*dt;
    double dt4 = dt3*dt;
    double q_pos = dt4/4*accel_var;
    double q_cross = dt3/2*accel_var;
    double q_vel = dt2*accel_var;
    for(int i=0;i<3;++i) {
        Q[i][i] = q_pos;
        Q[i][i+3] = q_cross;
        Q[i+3][i] = q_cross;
        Q[i+3][i+3] = q_vel;
    }
    return Q;
}

 Matrix6 directionalProcessNoise(double dt, double accel_var, const std::array<double,3>& vel) {
    Matrix6 Q = processNoise(dt, accel_var);
    double vx = vel[0];
    double vy = vel[1];
    double norm = std::sqrt(vx*vx + vy*vy);
    if(norm < 1e-6) return Q;
    double ux = vx / norm;
    double uy = vy / norm;

    Matrix6 T = identity6();
    // rotation for position
    T[0][0] = ux;  T[0][1] = uy;  T[0][2] = 0;
    T[1][0] = -uy; T[1][1] = ux;  T[1][2] = 0;
    T[2][2] = 1;
    // rotation for velocity
    T[3][3] = ux;  T[3][4] = uy;  T[3][5] = 0;
    T[4][3] = -uy; T[4][4] = ux;  T[4][5] = 0;
    T[5][5] = 1;

    Matrix6 Qloc = multiply6(multiply6(T, Q), transpose6(T));
    // increase cross-track noise (index 1 and 4)
    Qloc[1][1] *= 5.0;
    // directionalProcessNoise 함수 내에서
    Qloc[2][2] *= 2.0; // 고도 위치 노이즈 조정
    Qloc[5][5] *= 2.0; // 고도 속도 노이즈 조정
    Qloc[4][4] *= 5.0;
    // transform back
    Matrix6 Tt = transpose6(T);
    Q = multiply6(multiply6(Tt, Qloc), T);
    return Q;
}

// Kalman filter predict step
 void predict(State& state, Matrix6& P, double dt, double accel_var) {
    for(int i=0;i<3;++i) state.pos[i] += dt * state.vel[i];
    Matrix6 F = identity6();
    for(int i=0;i<3;++i) F[i][i+3]=dt;
    Matrix6 Ft = transpose6(F);
    Matrix6 temp = multiply6(F, P);
    P = add6(multiply6(temp, Ft), directionalProcessNoise(dt, accel_var, state.vel));
}

// Kalman filter update step
 void update(State& state, Matrix6& P, const std::array<double,3>& z, const Matrix3& R) {
    // Innovation
    std::array<double,3> y{};
    for(int i=0;i<3;++i) y[i] = z[i] - state.pos[i];

    // S = top-left 3x3 block of P plus R
    Matrix3 S{};
    for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
            S[i][j] = P[i][j] + R[i][j];
    Matrix3 Si = inverse3(S);

    // K = PHt * S^{-1}; PHt is first three columns of P
    Matrix6 K{}; // 6x3
    for(int i=0;i<6;++i)
        for(int j=0;j<3;++j) {
            double sum=0.0;
            for(int k=0;k<3;++k) sum += P[i][k] * Si[k][j];
            K[i][j] = sum;
        }

    // state = state + K*y
    for(int i=0;i<3;++i) {
        double inc = 0.0;
        for(int k=0;k<3;++k) inc += K[i][k]*y[k];
        state.pos[i] += inc;
    }
    for(int i=0;i<3;++i) {
        double inc = 0.0;
        for(int k=0;k<3;++k) inc += K[i+3][k]*y[k];
        state.vel[i] += inc;
    }

    // P = (I - K H) * P.  H selects first 3 state components.
    Matrix6 KH{};
    for(int i=0;i<6;++i)
        for(int j=0;j<6;++j)
            KH[i][j] = (j < 3) ? K[i][j] : 0.0;

    Matrix6 I = identity6();
    Matrix6 IKH = subtract6(I, KH);
    P = multiply6(IKH, P);
}

// RTS smoother
 void smooth(std::vector<StepData>& steps, double accel_var) {
    int n = steps.size();
    for(int k=n-2;k>=0;--k) {
        double dt = steps[k+1].dt;
        Matrix6 F = identity6();
        for(int i=0;i<3;++i) F[i][i+3]=dt;
        Matrix6 Ft = transpose6(F);
        Matrix6 P_pred = add6(multiply6(multiply6(F, steps[k].filtered_cov), Ft),
                                directionalProcessNoise(dt, accel_var, steps[k].filtered_state.vel));
        Matrix6 P_pred_inv = P_pred; // approximate inverse using Gauss-Jordan
        // We implement simple Gauss-Jordan inversion for 6x6
        const int N=6;
        Matrix6 A = P_pred;
        Matrix6 Inv = identity6();
        for(int i=0;i<N;i++) {
            double diag = A[i][i];
            if(std::abs(diag) < 1e-12) throw std::runtime_error("Singular P_pred");
            double invdiag = 1.0/diag;
            for(int j=0;j<N;j++) { A[i][j]*=invdiag; Inv[i][j]*=invdiag; }
            for(int k2=0;k2<N;k2++) if(k2!=i) {
                double f = A[k2][i];
                for(int j=0;j<N;j++) {
                    A[k2][j]-=f*A[i][j];
                    Inv[k2][j]-=f*Inv[i][j];
                }
            }
        }
        P_pred_inv = Inv;
        Matrix6 C = multiply6(multiply6(steps[k].filtered_cov, Ft), P_pred_inv);
        // smoothed state
        std::array<double,6> x_filt{steps[k].filtered_state.pos[0],steps[k].filtered_state.pos[1],steps[k].filtered_state.pos[2],steps[k].filtered_state.vel[0],steps[k].filtered_state.vel[1],steps[k].filtered_state.vel[2]};
        std::array<double,6> x_pred{steps[k+1].predicted_state.pos[0],steps[k+1].predicted_state.pos[1],steps[k+1].predicted_state.pos[2],steps[k+1].predicted_state.vel[0],steps[k+1].predicted_state.vel[1],steps[k+1].predicted_state.vel[2]};
        std::array<double,6> x_smooth_next{steps[k+1].filtered_state.pos[0],steps[k+1].filtered_state.pos[1],steps[k+1].filtered_state.pos[2],steps[k+1].filtered_state.vel[0],steps[k+1].filtered_state.vel[1],steps[k+1].filtered_state.vel[2]};
        // In smoothing, use x_{k+1|N}
        std::array<double,6> diff = subtract6(x_smooth_next, x_pred);
        std::array<double,6> corr = multiply6(C, diff);
        std::array<double,6> x_smooth = add6(x_filt, corr);
        for(int i=0;i<3;++i) {
            steps[k].filtered_state.pos[i] = x_smooth[i];
            steps[k].filtered_state.vel[i] = x_smooth[i+3];
        }
        Matrix6 P_smooth_next = steps[k+1].filtered_cov;
        Matrix6 temp = subtract6(P_smooth_next, P_pred);
        Matrix6 P_corr = multiply6(multiply6(C, temp), transpose6(C));
        steps[k].filtered_cov = add6(steps[k].filtered_cov, P_corr);
    }
}

#endif
