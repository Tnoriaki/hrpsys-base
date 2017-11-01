/* -*- coding:utf-8-unix; mode:c++; -*- */
#include "FootGuidedController.h"

using namespace hrp;
using namespace rats;

void foot_guided_control_base::set_param(const double _dz, const size_t N_max)
{
    // set param
    dz = _dz;
    xi = std::sqrt(g / dz);
    h = 1 + xi * dt;
    h_ = 1 - xi * dt;
    A <<
        1.0, dt,
        xi * xi * dt, 1.0;
    b <<
        0.0,
        -xi * xi * dt;
    Phi <<
        1.0, 1.0 / xi,
        1.0, -1.0 / xi;
    Phi_inv <<
        1.0, 1.0,
        xi, -xi;
    Phi_inv = 2.0 * Phi_inv;
    // calc h gain
    calc_refzmp_gain_list(N_max);
}

double foot_guided_control_base::calc_hgain(const std::size_t N)
{
    if (N > 0) {
        return (1 + h) * std::pow(h, N-1) / (1 - std::pow(h, 2*N));
    } else {
        return 0;
    }
}

void foot_guided_control_base::calc_hgain_list(const std::size_t N)
{
    hgain_list.resize(N);
    for (size_t i = 0; i < N; i++) {
        hgain_list(i) = calc_hgain(i);
    }
}

void foot_guided_control_base::calc_refzmp_gain_list(const std::size_t N)
{
    calc_hgain_list(N);
    refzmp_gain_list.resize(N);
    for (size_t i = 0; i < N; i++) {
        refzmp_gain_list(i) = hgain_list(i) * std::pow(h, N);
    }
}

void foot_guided_control_base::calc_u(const double remain_t, const double ref_dcm, const double ref_zmp)
{
    int N_remain = round(remain_t / dt);
    double hn = std::pow(h, N_remain);
    double hgain =  (hn / h) * (1 + h) / ( 1 - hn * hn );
    w_k = Phi * x_k;
    u_k = ref_zmp - hgain * (hn * (w_k(0) - ref_zmp) - (ref_dcm - ref_zmp));
}

void foot_guided_control_base::calc_u(const double remain_t, const double ref_dcm, const double ref_ccm, const double ref_zmp)
{
    // set param
    int N_remain = round(remain_t / dt);
    if (N_remain == 1) {
        calc_u(remain_t, ref_dcm, ref_zmp);
        return;
    }
    double ref_cm[2] = {ref_dcm, ref_ccm};
    double hn[2] = { std::pow(h, N_remain), std::pow(h_, N_remain) };
    double H[3] = {h * h, h_ * h_, h * h_};
    for (size_t i = 0; i < 3; i ++) H[i] = (1 - std::pow(H[i], N_remain)) / (1 - H[i]);
    double det = H[0] * H[1] - H[2] * H[2];
    double Hcoef[2];
    Hcoef[0] = std::pow(h, N_remain-1) * H[1] - std::pow(h_, N_remain-1) * H[2];
    Hcoef[1] = std::pow(h, N_remain-1) * H[2] - std::pow(h_, N_remain-1) * H[0];
    // calc u
    w_k = Phi * x_k;
    u_k = ref_zmp;
    for (size_t i = 0; i < 2; i++)
        u_k += 1 / (det * xi * dt) * Hcoef[i] * ( hn[i] * w_k(i) - ref_cm[i] + (1 - hn[i]) * ref_zmp );
}

// void foot_guided_control_base::calc_u_from_refzmp_list(const double remain_t, const double ref_dcm, const hrp::dvector ref_zmp_list )
// {
//     int N_remain = round(remain_t / dt);
//     w_k = Phi * x_k;
//     double sum = 0;
//     for (size_t i = 0; i < N_remain; i++) {
//         sum += xi * dt * std::pow(h, N_remain - i - 1) * ref_zmp_list(i);
//     }
//     u_k = ref_zmp_list(0) - hgain_list(N_remain) * (std::pow(h, N_remain) * w_k(0) - ref_dcm - sum);
// }

// assumed after calc_u
void foot_guided_control_base::calc_x_k()
{
    x_k = A * x_k + b * u_k;
}

void foot_guided_control_base::update_x_k(const double remain_t, const double ref_dcm, const double ref_zmp)
{
    calc_u(remain_t, ref_dcm, ref_zmp);
    calc_x_k();
}

void foot_guided_control_base::update_x_k(const double remain_t, const double ref_dcm, const double ref_ccm, const double ref_zmp)
{
    calc_u(remain_t, ref_dcm, ref_ccm, ref_zmp);
    calc_x_k();
}
